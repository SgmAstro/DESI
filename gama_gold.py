import os
import argparse
import runtime
import numpy           as np
import astropy.io.fits as fits

from   astropy.table   import Table
from   cosmo           import cosmo, distmod
from   gama_limits     import gama_field
from   cartesian       import cartesian, rotate


root   = os.environ['TILING_CATDIR']
fpath  = root + '/TilingCatv46.fits'
opath  = os.environ['GOLD_DIR'] + '/gama_gold.fits'

parser = argparse.ArgumentParser(description='Gen kE cat.')
parser.add_argument('--nooverwrite',  help='Do not overwrite outputs if on disk', action='store_true')

args   = parser.parse_args()

if args.nooverwrite:
    if os.path.isfile(opath):
        print('{} found on disk and overwrite forbidden (--nooverwrite).'.format(opath))
        exit(0)

dat    = Table.read(fpath)
dat    = Table(dat, masked=False)

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
z_cut      = (dat['ZGAMA'] > 0.039) & (dat['ZGAMA'] < 0.263)
r_cut      = (dat['R_PETRO'] > 12)
nq_cut     = (dat['NQ'] >= 3)

print(np.mean(sclass_cut))
print(np.mean(z_cut))
print(np.mean(r_cut))
print(np.mean(nq_cut))

dat = dat[sclass_cut & z_cut & r_cut & nq_cut]

dat['LUMDIST'] = cosmo.luminosity_distance(dat['ZGAMA'].data)
dat['DISTMOD'] = distmod(dat['ZGAMA'].data)

dat['FIELD']   = gama_field(dat['RA'], dat['DEC'])

xyz = cartesian(dat['RA'], dat['DEC'], dat['ZGAMA'])

dat['CARTESIAN_X'] = xyz[:,0]
dat['CARTESIAN_Y'] = xyz[:,1]
dat['CARTESIAN_Z'] = xyz[:,2]

xyz = rotate(dat['RA'], dat['DEC'], xyz)

dat['ROTCARTESIAN_X'] = xyz[:,0]
dat['ROTCARTESIAN_Y'] = xyz[:,1]
dat['ROTCARTESIAN_Z'] = xyz[:,2]

# Randomise rows.
idx = np.arange(len(dat))
idx = np.random.choice(idx, size=len(idx), replace=False)

dat = dat[idx]

print('Writing {}.'.format(os.environ['GOLD_DIR'] + '/gama_gold.fits'))

if not os.path.isdir(os.environ['GOLD_DIR']):
    print('Creating {}'.format(os.environ['GOLD_DIR']))

    os.makedirs(os.environ['GOLD_DIR'])

# 113687 vs TMR 80922.
dat.meta['GOLD_NGAL'] = len(dat)
dat.pprint()
dat.write(opath, format='fits', overwrite=True)

idx   = np.random.choice(np.arange(len(dat)), 5000, replace=False)
dat   = dat[idx]

dat.write(os.environ['CODE_ROOT'] + '/data/gama_gold_dryrun.fits', format='fits', overwrite=True)
