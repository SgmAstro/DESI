import os
import argparse
import fitsio

from   astropy.table import Table
from   ddp           import get_ddps, tmr_DDP1, tmr_DDP2, tmr_DDP3


parser = argparse.ArgumentParser(description='Gen ddp cat.')
parser.add_argument('-d', '--dryrun', help='Dryrun.', action='store_true')

args   = parser.parse_args()
dryrun = args.dryrun

fpath  = os.environ['GOLD_DIR'] + '/gama_gold_zmax.fits'

if dryrun:
    fpath = fpath.replace('.fits', '_dryrun.fits')

print('Reading: {}'.format(fpath))
    
dat    = Table.read(fpath)
Area   = dat.meta['AREA']

print('Retrieved Area: {}'.format(Area))
print('Judging DDP.')

dat['DDP'], zlims = get_ddps(Area, dat['MALL_0P0'], dat['ZGAMA'])
dat.meta.update(zlims)
dat.meta.update({'TMR_DDP1': str(tmr_DDP1),\
                {'TMR_DDP2': str(tmr_DDP2),\
                {'TMR_DDP3': str(tmr_DDP1)})

print(zlims)

opath = fpath.replace('zmax', 'ddp')

print('Writing: {}'.format(opath))

dat.write(opath, format='fits', overwrite=True)
