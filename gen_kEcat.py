import os

from   astropy.table import Table
from smith_kcorr import GAMA_KCorrection


root = os.environ['CSCRATCH'] + '/norberg/'
fpath = root + '/GAMA4/gama_gold.fits'

dat = Table.read(fpath)
dat.pprint()



#
med_smith_color = 0.603

kcorr_r = GAMA_KCorrection(band='R')
kcorr_g = GAMA_KCorrection(band='G')

# dat['RESTGMR_0P1'] = med_smith_color * np.ones(len(dat))

dat['RKCORR_0P1']  = GAMA_KCorrection(band='R').k(dat['ZGAMA'], )
