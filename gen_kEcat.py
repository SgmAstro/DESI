import os
import argparse
import runtime
import numpy as np

from   astropy.table import Table
from   smith_kcorr   import GAMA_KCorrection
from   rest_gmr      import smith_rest_gmr
from   tmr_ecorr     import tmr_ecorr, tmr_q
from   abs_mag       import abs_mag

from   gama_limits   import gama_field, gama_fields
from   desi_fields   import desi_fields
from   findfile      import findfile, fetch_fields, overwrite_check


np.random.seed(314)

parser  = argparse.ArgumentParser(description='Gen kE cat.')
parser.add_argument('-d', '--dryrun', help='Dryrun.', action='store_true')
parser.add_argument('-s', '--survey', help='Select survey', default='gama')
parser.add_argument('--nooverwrite',  help='Do not overwrite outputs if on disk', action='store_true')

args        = parser.parse_args()
dryrun      = args.dryrun
survey      = args.survey.lower()
nooverwrite = args.nooverwrite

zsurv   = f'z{survey}'.upper()

root    = os.environ['GOLD_DIR']

fpath   = findfile(ftype='gold', dryrun=False, survey=survey)
opath   = findfile(ftype='kE',   dryrun=dryrun, survey=survey)

if args.nooverwrite:
  overwrite_check(opath)

fields  = fetch_fields(survey)
    
dat     = Table.read(fpath)
dat.pprint()

rest_gmr_0p1, rest_gmr_0p1_warn = smith_rest_gmr(dat[zsurv], dat['GMR'])

dat['REST_GMR_0P1']      = rest_gmr_0p1
dat['REST_GMR_0P1_WARN'] = rest_gmr_0p1_warn.astype(np.int32)

kcorr_r = GAMA_KCorrection(band='R')
kcorr_g = GAMA_KCorrection(band='G')

dat['REST_GMR_0P1_INDEX'] = kcorr_r.rest_gmr_index(dat['REST_GMR_0P1'], kcoeff=False)

dat['KCORR_R0P1'] = kcorr_r.k(dat[zsurv], dat['REST_GMR_0P1'])
dat['KCORR_G0P1'] = kcorr_g.k(dat[zsurv], dat['REST_GMR_0P1'])

dat['KCORR_R0P0'] = kcorr_r.k_nonnative_zref(0.0, dat[zsurv], dat['REST_GMR_0P1'])
dat['KCORR_G0P0'] = kcorr_g.k_nonnative_zref(0.0, dat[zsurv], dat['REST_GMR_0P1'])

dat['REST_GMR_0P0'] = dat['GMR'] - (dat['KCORR_G0P0'] - dat['KCORR_R0P0'])

dat['Q_COLOR_0P0'] = tmr_q(dat['REST_GMR_0P0'], aall=False)

dat['EQ_ALL_0P0']   = tmr_ecorr(dat[zsurv], dat['REST_GMR_0P0'], aall=True)
dat['EQ_COLOR_0P0'] = tmr_ecorr(dat[zsurv], dat['REST_GMR_0P0'], aall=False)

dat['MALL_0P0']     = abs_mag(dat['DETMAG'], dat['DISTMOD'], dat['KCORR_R0P0'], dat['EQ_ALL_0P0'])
dat['MCOLOR_0P0']   = abs_mag(dat['DETMAG'], dat['DISTMOD'], dat['KCORR_R0P0'], dat['EQ_COLOR_0P0'])
dat['MQZERO_0P0']   = abs_mag(dat['DETMAG'], dat['DISTMOD'], dat['KCORR_R0P0'], np.zeros_like(dat['EQ_ALL_0P0']))

dat['Z_THETA_QALL']   = dat['DISTMOD'] + dat['KCORR_R0P0'] + dat['EQ_ALL_0P0']
dat['Z_THETA_QZERO']  = dat['DISTMOD'] + dat['KCORR_R0P0'] + np.zeros_like(dat['EQ_ALL_0P0'])
dat['Z_THETA_QCOLOR'] = dat['DISTMOD'] + dat['KCORR_R0P0'] + dat['EQ_COLOR_0P0']

##  ----  DDP  ----
##  Note:  assumes median rest-frame colour and QALL.
dat['DDPKCORR_R0P1'] = kcorr_r.k(dat[zsurv], dat['REST_GMR_0P1'], median=True)
dat['DDPKCORR_G0P1'] = kcorr_g.k(dat[zsurv], dat['REST_GMR_0P1'], median=True)

dat['DDPKCORR_R0P0'] = kcorr_r.k_nonnative_zref(0.0, dat[zsurv], dat['REST_GMR_0P1'], median=True)
dat['DDPKCORR_G0P0'] = kcorr_g.k_nonnative_zref(0.0, dat[zsurv], dat['REST_GMR_0P1'], median=True)

dat['DDPMALL_0P0']   = abs_mag(dat['DETMAG'], dat['DISTMOD'], dat['DDPKCORR_R0P0'], dat['EQ_ALL_0P0'])

dat.meta['IMMUTABLE'] = 'False'

dat.pprint()

print('Writing {}.'.format(opath))

dat.write(opath, format='fits', overwrite=True)
