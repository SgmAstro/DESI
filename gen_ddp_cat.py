import os
import sys
import argparse
import runtime
import fitsio
import numpy         as     np

from   astropy.table import Table
from   ddp           import get_ddps, tmr_DDP1, tmr_DDP2, tmr_DDP3, _initialise_ddplimits
from   ddp_limits    import limiting_curve_path
from   findfile      import findfile, overwrite_check, write_desitable
from   bitmask       import lumfn_mask, consv_mask, update_bit
from   config        import Configuration


parser = argparse.ArgumentParser(description='Gen ddp cat.')
parser.add_argument('--log', help='Create a log file of stdout.', action='store_true')
parser.add_argument('-d', '--dryrun', help='Dryrun.', action='store_true')
parser.add_argument('-s', '--survey', help='Select survey', default='gama')
parser.add_argument('--config',       help='Path to configuration file', type=str, default=findfile('config'))
parser.add_argument('--nooverwrite',  help='Do not overwrite outputs if on disk', action='store_true')

args   = parser.parse_args()
log    = args.log
dryrun = args.dryrun
survey = args.survey

config = Configuration(args.config)
config.update_attributes('ddp', args)
config.write()

fpath  = findfile(ftype='zmax', dryrun=dryrun, survey=survey)
opath  = findfile(ftype='ddp',  dryrun=dryrun, survey=survey)

if log:
    logfile = findfile(ftype='ddp', dryrun=False, survey=survey, log=True)

    print(f'Logging to {logfile}')
    
    sys.stdout = open(logfile, 'w')

if args.nooverwrite:
    overwrite_check(opath)

print('Reading: {}'.format(fpath))
    
dat    = Table.read(fpath)
Area   = dat.meta['AREA']

print('Retrieved Area: {}'.format(Area))
print('Judging DDP.')

dat['DDP'], dat['DDPZLIMS'], zlims = get_ddps(Area, dat['DDPMALL_0P0'], dat['ZSURV'], survey)

dat['STEPWISE_FAINTLIM_0P0']  = -99.
dat['STEPWISE_BRIGHTLIM_0P0'] = -99.

for color_idx in np.unique(dat['REST_GMR_0P1_INDEX']):
    isin                                 = (dat['REST_GMR_0P1_INDEX'] == color_idx)

    # Stepwise limiting magnitudes for a given zref=0.1 color and redshift.
    fidx, faint_limit_path               = limiting_curve_path(survey, dat.meta['RLIM'], 'QCOLOR', gmr_0P1_idx=color_idx, gmr_0P1=None, gmr_0P0=None, debug=False)
    bidx, bright_limit_path              = limiting_curve_path(survey, dat.meta['RMAX'], 'QCOLOR', gmr_0P1_idx=color_idx, gmr_0P1=None, gmr_0P0=None, debug=False)

    _, bright_curve_r, _, faint_curve_r  = _initialise_ddplimits(bidx, fidx, survey=survey, Mcol='M0P0_QCOLOR')

    dat['STEPWISE_FAINTLIM_0P0'][isin]   =  faint_curve_r(dat['ZSURV'][isin]) 
    dat['STEPWISE_BRIGHTLIM_0P0'][isin]  = bright_curve_r(dat['ZSURV'][isin])

update_bit(dat['IN_D8LUMFN'], lumfn_mask, 'DDP1ZLIM', dat['DDPZLIMS'][:,0] == 0)

dat.meta.update(zlims)
dat.meta.update({'TMR_DDP1': str(tmr_DDP1),\
                 'TMR_DDP2': str(tmr_DDP2),\
                 'TMR_DDP3': str(tmr_DDP3)})

print(zlims)

print('Writing: {}'.format(opath))

write_desitable(opath, dat)

if log:
    sys.stdout.close()
