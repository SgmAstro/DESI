import os
import numpy as np

from astropy.table import Table
from smith_kcorr import GAMA_KCorrection
from rest_gmr import smith_rest_gmr
from tmr_ecorr import tmr_ecorr, tmr_q
from abs_mag import abs_mag

#ngal=5000
ngal=100000
nproc=4

root = os.environ['CSCRATCH'] + '/norberg/'
fpath = root + '/GAMA4/gama_gold.fits'

dat = Table.read(fpath)
dat.pprint()

dat = Table(np.random.choice(dat, ngal))

dat['GMR'] = dat['GMAG_DRED_SDSS'] - dat['RMAG_DRED_SDSS']

rest_gmr_0p1, rest_gmr_0p1_warn = smith_rest_gmr(dat['ZGAMA'], dat['GMR'])

dat['REST_GMR_0P1'] = rest_gmr_0p1
dat['REST_GMR_0P1_WARN'] = rest_gmr_0p1_warn.astype(np.int)

kcorr_r = GAMA_KCorrection(band='R')
kcorr_g = GAMA_KCorrection(band='G')

dat['REST_GMR_0P1_INDEX'] = kcorr_r.rest_gmr_index(dat['REST_GMR_0P1'], kcoeff=False)

dat['KCORR_R0P1'] = kcorr_r.k(dat['ZGAMA'], dat['REST_GMR_0P1'])
dat['KCORR_G0P1'] = kcorr_g.k(dat['ZGAMA'], dat['REST_GMR_0P1'])

dat['KCORR_R0P0'] = kcorr_r.k_nonnative_zref(0.0, dat['ZGAMA'], dat['REST_GMR_0P1'])
dat['KCORR_G0P0'] = kcorr_g.k_nonnative_zref(0.0, dat['ZGAMA'], dat['REST_GMR_0P1'])

dat['REST_GMR_0P0'] = dat['GMR'] - (dat['KCORR_G0P0'] - dat['KCORR_R0P0'])

dat['Q_COLOR_0P0'] = tmr_q(dat['ZGAMA'], dat['REST_GMR_0P0'], aall=False)

dat['EQ_ALL_0P0']   = tmr_ecorr(dat['ZGAMA'], dat['REST_GMR_0P0'], aall=True)
dat['EQ_COLOR_0P0']   = tmr_ecorr(dat['ZGAMA'], dat['REST_GMR_0P0'], aall=False)

dat['MALL_0P0'] = abs_mag(dat['R_PETRO'], dat['DISTMOD'], dat['KCORR_R0P0'], dat['EQ_ALL_0P0'])
dat['MCOLOR_0P0'] = abs_mag(dat['R_PETRO'], dat['DISTMOD'], dat['KCORR_R0P0'], dat['EQ_COLOR_0P0'])

dat['Z_THETA_QALL'] = dat['DISTMOD'] + dat['KCORR_R0P0'] + dat['EQ_ALL_0P0']
dat['Z_THETA_QCOLOR'] = dat['DISTMOD'] + dat['KCORR_R0P0'] + dat['EQ_COLOR_0P0']

dat.pprint()

dat.write(root + '/GAMA4/gama_gold_kE_{:d}k.fits'.format(np.int(ngal / 1000.)), format='fits', overwrite=True)
