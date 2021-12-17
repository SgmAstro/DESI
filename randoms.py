import os
import numpy as np

from   cosmo import cosmo, volcom
from   scipy.interpolate import interp1d
from   gama_limits import gama_limits
from   astropy.table import Table
from   cartesian import cartesian


np.random.seed(314)

realz = 0
field = 'G9'
Area = 60. 

dz   = 1.e-4

zmin = 0.0
zmax = 0.3

Vmin = volcom(zmin, Area) 
Vmax = volcom(zmax, Area)

# Assumse one gama field, of 60. sq. deg. 
vol  = Vmax - Vmin
rand_density = 5.e-1

nrand = np.int(np.ceil(vol * rand_density))

print(vol, rand_density, nrand / 1.e6)

boundary_percent = 1.

opath = os.environ['CSCRATCH'] + '/desi/BGS/Sam/randoms_{}_{:d}.fits'.format(field, realz)

zs   = np.arange(0.0, zmax+dz, dz)
Vs   = volcom(zs, Area) 

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

print('Applying rotation.')

xyz = cartesian(ras, decs, zs)

ras = ras.astype(np.float32)
decs= decs.astype(np.float32)
zs = zs.astype(np.float32)
Vdraws = Vdraws.astype(np.float32)
xyz=xyz.astype(np.float32)

randoms = Table(np.c_[ras, decs, zs, Vdraws], names=['RANDOM_RA', 'RANDOM_DEC', 'Z', 'V'])
randoms['RANDID'] = np.arange(len(randoms))
randoms['FIELD'] = field

randoms['CARTESIAN_X'] = xyz[:,0]
randoms['CARTESIAN_Y'] = xyz[:,1]
randoms['CARTESIAN_Z'] = xyz[:,2]

print('Applying boundary.')

randoms['IS_BOUNDARY'] = 0

randoms['IS_BOUNDARY'][randoms['RANDOM_RA']  > np.percentile(randoms['RANDOM_RA'], 100. - boundary_percent)] = 1
randoms['IS_BOUNDARY'][randoms['RANDOM_RA']  < np.percentile(randoms['RANDOM_RA'],  boundary_percent)] = 1

randoms['IS_BOUNDARY'][randoms['RANDOM_DEC'] > np.percentile(randoms['RANDOM_DEC'], 100. - boundary_percent)] = 1
randoms['IS_BOUNDARY'][randoms['RANDOM_DEC'] < np.percentile(randoms['RANDOM_DEC'], boundary_percent)] = 1

randoms['IS_BOUNDARY'][randoms['V'] >= np.percentile(randoms['V'], 100. - boundary_percent)] = 1
randoms['IS_BOUNDARY'][randoms['V'] <= np.percentile(randoms['V'],  boundary_percent)] = 1

randoms.meta = {'ZMIN': zmin, 'ZMAX': zmax, 'DZ': dz, 'NRAND': nrand, 'FIELD': field, 'Area': Area, 'BOUND_PERCENT': boundary_percent,\
                'VOL': vol, 'RAND_DENS': rand_density, 'VOL8': (4./3.)*np.pi*(8.**3.)}

randoms.meta['NRAND8'] = randoms.meta['VOL8'] * randoms.meta['RAND_DENS']
randoms.meta['NRAND8_PERR'] =np.sqrt(randoms.meta['NRAND8'])

print(randoms.meta)

print('Writing.')

# randoms.write(opath, format='fits', overwrite=True)
