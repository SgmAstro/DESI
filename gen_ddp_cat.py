import os
import argparse
import runtime
import fitsio

from   astropy.table import Table
from   ddp           import get_ddps, tmr_DDP1, tmr_DDP2, tmr_DDP3
from   findfile      import findfile, overwrite_check
from   bitmask       import lumfn_mask, consv_mask


parser = argparse.ArgumentParser(description='Gen ddp cat.')
parser.add_argument('--log', help='Create a log file of stdout.', action='store_true')
parser.add_argument('-d', '--dryrun', help='Dryrun.', action='store_true')
parser.add_argument('-s', '--survey', help='Select survey', default='gama')
parser.add_argument('--nooverwrite',  help='Do not overwrite outputs if on disk', action='store_true')

args   = parser.parse_args()
dryrun = args.dryrun
survey = args.survey

fpath  = findfile(ftype='zmax', dryrun=dryrun, survey=survey)
opath  = findfile(ftype='ddp',  dryrun=dryrun, survey=survey)

if args.nooverwrite:
    overwrite_check(opath)

print('Reading: {}'.format(fpath))
    
dat    = Table.read(fpath)
Area   = dat.meta['AREA']

print('Retrieved Area: {}'.format(Area))
print('Judging DDP.')

dat['DDP'], dat['DDPZLIMS'], zlims, _ = get_ddps(Area, dat['DDPMALL_0P0'], dat['ZSURV'], survey)

dat['IN_D8LUMFN'] += (dat['DDPZLIMS'][:,0] == 0) * lumfn_mask.DDP1ZLIM

dat.meta.update(zlims)
dat.meta.update({'TMR_DDP1': str(tmr_DDP1),\
                 'TMR_DDP2': str(tmr_DDP2),\
                 'TMR_DDP3': str(tmr_DDP3)})

print(zlims)

print('Writing: {}'.format(opath))

dat.write(opath, format='fits', overwrite=True)
