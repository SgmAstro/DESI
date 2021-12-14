import os

from astropy.table import Table
from smith_kcorr import GAMA_KCorrection
from rest_gmr import smith_rest_gmr
from tmr_ecorr import tmr_ecorr, tmr_q


root = os.environ['CSCRATCH'] + '/norberg/'
fpath = root + '/GAMA4/gama_gold.fits'

dat = Table.read(fpath)
dat.pprint()

dat = dat[:5000]

dat['GMR'] = dat['GMAG_DRED_SDSS'] - dat['RMAG_DRED_SDSS']

rest_gmr_0p1, rest_gmr_0p1_warn = smith_rest_gmr(dat['ZGAMA'], dat['GMR'])

dat['REST_GMR_0P1'] = rest_gmr_0p1
dat['REST_GMR_0P1_WARN'] = rest_gmr_0p1_warn

kcorr_r = GAMA_KCorrection(band='R')
kcorr_g = GAMA_KCorrection(band='G')

dat['KCORR_R0P1'] = kcorr_r.k(dat['ZGAMA'], dat['REST_GMR_0P1'])
dat['KCORR_G0P1'] = kcorr_g.k(dat['ZGAMA'], dat['REST_GMR_0P1'])

dat['KCORR_R0P0'] = kcorr_r.k_nonnative_zref(0.0, dat['ZGAMA'], dat['REST_GMR_0P1'])
dat['KCORR_G0P0'] = kcorr_g.k_nonnative_zref(0.0, dat['ZGAMA'], dat['REST_GMR_0P1'])

dat['REST_GMR_0P0'] = dat['GMR'] - (dat['KCORR_G0P0'] - dat['KCORR_R0P0'])

dat['Q_COLOR_0P0'] = tmr_q(dat['ZGAMA'], dat['REST_GMR_0P0'], all=False)

dat['EQ_ALL_0P0']   = tmr_ecorr(dat['ZGAMA'], dat['REST_GMR_0P0'], all=True)
dat['EQ_COLOR_0P0']   = tmr_ecorr(dat['ZGAMA'], dat['REST_GMR_0P0'], all=False)

dat.pprint()
