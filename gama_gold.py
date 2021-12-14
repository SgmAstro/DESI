import os
import numpy as np
import astropy.io.fits as fits

from astropy.table import Table
from cosmo import cosmo


root = os.environ['CSCRATCH'] + '/norberg/'
fpath = root + '/TilingCatv46.fits'

dat   = Table.read(fpath)

# print(dat.dtype.names)

dat['ZGAMA'] = dat['Z']

del dat['Z']

for band in 'UGRIZ':
    dat['{}MAG_DRED_SDSS'.format(band)] = dat['{}_MODEL'.format(band)]

    del dat['{}_MODEL'.format(band)]
    
minimal_cols = ['CATAID', 'OBJID', 'RA', 'DEC', 'R_PETRO', 'ZGAMA', 'NQ', 'SPECID', 'SURVEY_CLASS']

for band in ['U', 'G', 'R', 'I', 'Z']:
    minimal_cols += ['{}MAG_DRED_SDSS'.format(band)]

# Minimal catalogue.
dat = dat[minimal_cols]
dat.pprint()

# 'SURVEY_CLASS' < 4 for GAMA-II (rpet <= 19.8 by extension.
# 0.039 < z < 0.263, DDP1 sample.	
# r=12 bright cut;
# 1 cat. per field (G9, 12, 15).

sclass_cut = dat['SURVEY_CLASS'] >= 4
z_cut = (dat['ZGAMA'] > 0.039) & (dat['ZGAMA'] < 0.263)
r_cut = (dat['R_PETRO'] > 12)
nq_cut = dat['NQ'] >= 3

print(np.mean(sclass_cut))
print(np.mean(z_cut))
print(np.mean(r_cut))
print(np.mean(nq_cut))

dat = dat[sclass_cut & z_cut & r_cut & nq_cut]
dat.pprint()

dat['LUMDIST'] = cosmo.luminosity_distance(dat['ZGAMA'])
dat['DISTMOD'] = 5. * np.log10(1.e6 * dat['LUMDIST'] / 10.)

print('Writing {}'.format(root + '/GAMA4/gama_gold.fits'))

dat.write(root + '/GAMA4/gama_gold.fits', format='fits', overwrite=True)
