import os
import argparse
import numpy as np

from astropy.table import Table
from smith_kcorr import GAMA_KCorrection
from rest_gmr import smith_rest_gmr
from tmr_ecorr import tmr_ecorr, tmr_q
from abs_mag import abs_mag


np.random.seed(314)

parser = argparse.ArgumentParser(description='Gen kE cat.')
parser.add_argument('-d', '--dryrun', help='Dryrun.', action='store_true')

args = parser.parse_args()
dryrun = args.dryrun

root = os.environ['GOLD_DIR']
fpath = root + '/gama_gold.fits'

dat = Table.read(fpath)
dat.pprint()

opath=root + '/gama_gold_kE.fits'

if dryrun:
  idx   = np.random.choice(np.arange(len(dat)), 500, replace=False)
  
  dat   = dat[idx]
  opath = opath.replace('.fits', '_dryrun.fits')                                                     

  
dat['GMR'] = dat['GMAG_DRED_SDSS'] - dat['RMAG_DRED_SDSS']

rest_gmr_0p1, rest_gmr_0p1_warn = smith_rest_gmr(dat['ZGAMA'], dat['GMR'])

dat['REST_GMR_0P1'] = rest_gmr_0p1
dat['REST_GMR_0P1_WARN'] = rest_gmr_0p1_warn.astype(np.int32)

kcorr_r = GAMA_KCorrection(band='R')
kcorr_g = GAMA_KCorrection(band='G')

dat['REST_GMR_0P1_INDEX'] = kcorr_r.rest_gmr_index(dat['REST_GMR_0P1'], kcoeff=False)

dat['KCORR_R0P1'] = kcorr_r.k(dat['ZGAMA'], dat['REST_GMR_0P1'])
dat['KCORR_G0P1'] = kcorr_g.k(dat['ZGAMA'], dat['REST_GMR_0P1'])

dat['KCORR_R0P0'] = kcorr_r.k_nonnative_zref(0.0, dat['ZGAMA'], dat['REST_GMR_0P1'])
dat['KCORR_G0P0'] = kcorr_g.k_nonnative_zref(0.0, dat['ZGAMA'], dat['REST_GMR_0P1'])

dat['REST_GMR_0P0'] = dat['GMR'] - (dat['KCORR_G0P0'] - dat['KCORR_R0P0'])

dat['Q_COLOR_0P0'] = tmr_q(dat['REST_GMR_0P0'], aall=False)

dat['EQ_ALL_0P0']   = tmr_ecorr(dat['ZGAMA'], dat['REST_GMR_0P0'], aall=True)
dat['EQ_COLOR_0P0']   = tmr_ecorr(dat['ZGAMA'], dat['REST_GMR_0P0'], aall=False)

dat['MALL_0P0'] = abs_mag(dat['R_PETRO'], dat['DISTMOD'], dat['KCORR_R0P0'], dat['EQ_ALL_0P0'])
dat['MCOLOR_0P0'] = abs_mag(dat['R_PETRO'], dat['DISTMOD'], dat['KCORR_R0P0'], dat['EQ_COLOR_0P0'])

dat['Z_THETA_QALL'] = dat['DISTMOD'] + dat['KCORR_R0P0'] + dat['EQ_ALL_0P0']
dat['Z_THETA_QCOLOR'] = dat['DISTMOD'] + dat['KCORR_R0P0'] + dat['EQ_COLOR_0P0']

dat.pprint()

print('Writing {}.'.format(opath))

dat.write(opath, format='fits', overwrite=True)
