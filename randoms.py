import os
import numpy as np

from cosmo import cosmo
from scipy.interpolate import interp1d
from gama_limits import gama_limits
from astropy.table import Table
from cartesian import cartesian
from fsky import fsky


np.random.seed(314)

field = 'G9'

dz   = 1.e-4

zmin = 0.0
zmax = 0.6

# Assumse one gama field, of 60. sq. deg. 
vol  = fsky(60.) * (cosmo.comoving_volume(zmax).value - cosmo.comoving_volume(zmin).value)
rand_density = 5.e-2

nrand = np.int(np.ceil(vol * rand_density))

print(vol, rand_density, nrand / 1.e6)

boundary_percent = 0.5

fpath = os.environ['CSCRATCH'] + '/desi/BGS/Sam/randoms.fits'

Vmin = cosmo.comoving_volume(zmin).value
Vmax = cosmo.comoving_volume(zmax).value

zs   = np.arange(0.0, zmax+dz, dz)
Vs   = cosmo.comoving_volume(zs).value

# https://docs.scipy.org/doc/scipy/reference/generated/scipy.interpolate.interp1d.html
Vz   = interp1d(Vs, zs, kind='linear', copy=True, bounds_error=True, fill_value=np.NaN, assume_sorted=False)

ra_min  = gama_limits[field]['ra_min']
ra_max  = gama_limits[field]['ra_max']

dec_min = gama_limits[field]['dec_min']
dec_max = gama_limits[field]['dec_max']
    
Vdraws  = np.random.uniform(0., 1., nrand)
Vdraws  = Vmin + Vdraws * (Vmax - Vmin)

zs      = Vz(Vdraws)

ctheta_min = np.cos(np.pi/2. - np.radians(dec_min))
ctheta_max = np.cos(np.pi/2  - np.radians(dec_max))
    
cos_theta = np.random.uniform(ctheta_min, ctheta_max, nrand)
theta     = np.arccos(cos_theta)
decs      = np.pi/2. - theta
decs      = np.degrees(decs)

ras       = np.random.uniform(ra_min, ra_max, nrand)

print('Solved {:d} for field {}'.format(nrand, field))

randoms = Table(np.c_[ras, decs, zs, Vdraws], names=['RANDOM_RA', 'RANDOM_DEC', 'Z', 'V'])
randoms['RANDID'] = np.arange(len(randoms))
randoms['FIELD'] = field

xyz                    = cartesian(randoms['RANDOM_RA'].data, randoms['RANDOM_DEC'].data, randoms['Z'].data)
randoms['CARTESIAN_X'] = xyz[:,0]
randoms['CARTESIAN_Y'] = xyz[:,1]
randoms['CARTESIAN_Z'] = xyz[:,2]

randoms['IS_BOUNDARY'] = 0

randoms['IS_BOUNDARY'][randoms['RANDOM_RA']  > np.percentile(randoms['RANDOM_RA'], 100. - boundary_percent)] = 1
randoms['IS_BOUNDARY'][randoms['RANDOM_RA']  < np.percentile(randoms['RANDOM_RA'],  boundary_percent)] = 1

randoms['IS_BOUNDARY'][randoms['RANDOM_DEC'] > np.percentile(randoms['RANDOM_DEC'], 100. - boundary_percent)] = 1
randoms['IS_BOUNDARY'][randoms['RANDOM_DEC'] < np.percentile(randoms['RANDOM_DEC'], boundary_percent)] = 1

randoms['IS_BOUNDARY'][randoms['V'] >= np.percentile(randoms['V'], 100. - boundary_percent)] = 1
randoms['IS_BOUNDARY'][randoms['V'] <= np.percentile(randoms['V'],  boundary_percent)] = 1

randoms.meta = {'ZMIN': zmin, 'ZMAX': zmax, 'DZ': dz, 'NRAND': nrand}

randoms.write(fpath, format='fits', overwrite=True)
