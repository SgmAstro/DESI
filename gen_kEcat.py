import os

from   astropy.table import Table
from smith_kcorr import GAMA_KCorrection
from rest_gmr import smith_rest_gmr

root = os.environ['CSCRATCH'] + '/norberg/'
fpath = root + '/GAMA4/gama_gold.fits'

dat = Table.read(fpath)
dat.pprint()

dat['GMR'] = dat['GMAG_DRED_SDSS'] - dat['RMAG_DRED_SDSS']

dat['REST_GMR'] = smith_rest_gmr(dat['ZGAMA'], dat['GMR'])

'''
#
med_smith_color = 0.603

kcorr_r = GAMA_KCorrection(band='R')
kcorr_g = GAMA_KCorrection(band='G')

# dat['RESTGMR_0P1'] = med_smith_color * np.ones(len(dat))

dat['RKCORR_0P1']  = GAMA_KCorrection(band='R').k(dat['ZGAMA'], )
'''
