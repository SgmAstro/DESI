import os
import numpy as np
import astropy.io.fits as fits

from astropy.table import Table
from cosmo import cosmo, distmod
from gama_limits import gama_field
from cartesian import cartesian


root  = os.environ['TILING_CATDIR']
fpath = root + '/TilingCatv46.fits'

dat   = Table.read(fpath)
dat   = Table(dat, masked=False)

# TODO: Inherit this from somewhere sensible.
dat.meta['AREA'] = 180.

# print(dat.dtype.names)

dat.rename_column('Z', 'ZGAMA')

for band in 'UGRIZ':
    dat.rename_column('{}_MODEL'.format(band), '{}MAG_DRED_SDSS'.format(band))
    
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

sclass_cut = (dat['SURVEY_CLASS'] >= 4)
z_cut = (dat['ZGAMA'] > 0.039) & (dat['ZGAMA'] < 0.263)
r_cut = (dat['R_PETRO'] > 12)
nq_cut = (dat['NQ'] >= 3)

print(np.mean(sclass_cut))
print(np.mean(z_cut))
print(np.mean(r_cut))
print(np.mean(nq_cut))

dat = dat[sclass_cut & z_cut & r_cut & nq_cut]
dat.pprint()

dat['LUMDIST'] = cosmo.luminosity_distance(dat['ZGAMA'].data)
dat['DISTMOD'] = distmod(dat['ZGAMA'].data)

dat['FIELD']   = gama_field(dat['RA'], dat['DEC'])

xyz = cartesian(dat['RA'], dat['DEC'], dat['ZGAMA'])

dat['CARTESIAN_X'] = xyz[:,0]
dat['CARTESIAN_Y'] = xyz[:,1]
dat['CARTESIAN_Z'] = xyz[:,2]

# Randomise rows.
idx = np.arange(len(dat))
idx = np.random.choice(idx, size=len(idx), replace=False)

dat = dat[idx]

print('Writing {}.'.format(os.environ['GOLD_DIR'] + '/gama_gold.fits'))

if not os.path.isdir(os.environ['GOLD_DIR']):
    print('Creating {}'.format(os.environ['GOLD_DIR']))

    os.makedirs(os.environ['GOLD_DIR'])

dat.write(os.environ['GOLD_DIR'] + '/gama_gold.fits', format='fits', overwrite=True)
