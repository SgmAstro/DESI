import os
import time
import numpy as np
import argparse

from   cosmo import cosmo, volcom
from   scipy.interpolate import interp1d
from   gama_limits import gama_limits
from   astropy.table import Table
from   cartesian import cartesian
from   gama_limits import gama_field
from   runtime import calc_runtime


np.random.seed(314)

parser  = argparse.ArgumentParser(description='Select GAMA field.')
parser.add_argument('-f', '--field',  type=str, help='select equatorial GAMA field: G9, G12, G15', required=True)
parser.add_argument('-d', '--dryrun', help='Dryrun.', action='store_true')
parser.add_argument('--prefix',       help='filename prefix', default='randoms')
parser.add_argument('--nooverwrite',  help='Do not overwrite outputs if on disk', action='store_true')

# Defaults to GAMA Gold limits. 
parser.add_argument('--zmin', type=np.float32, help='Minimum redshift limit', default=0.039)
parser.add_argument('--zmax', type=np.float32, help='Maximum redshift limit', default=0.263)

args    = parser.parse_args()
field   = args.field.upper()
dryrun  = args.dryrun
zmin    = args.zmin
zmax    = args.zmax
prefix  = args.prefix 


start   = time.time()

realz   = 0
Area    = 60. 

dz      = 1.e-4

Vmin    = volcom(zmin, Area) 
Vmax    = volcom(zmax, Area)

# Assumse one gama field, of 60. sq. deg. 
vol          = Vmax - Vmin
rand_density = 1.

nrand        = np.int64(np.ceil(vol * rand_density))

opath        = os.environ['RANDOMS_DIR'] + '/{}_{}_{:d}.fits'.format(prefix, field, realz)

if dryrun:
    nrand = 500
    opath = opath.replace('.fits', '_dryrun.fits')

if args.nooverwrite:
    if os.path.isfile(opath):
        print('{} found on disk and overwrite forbidden (--nooverwrite).'.format(opath))
        exit(0)

if not os.path.isdir(os.environ['RANDOMS_DIR']):
    print('Creating {}'.format(os.environ['RANDOMS_DIR']))

    os.makedirs(os.environ['RANDOMS_DIR'])

    
print('Volume [1e6]: {:.2f}; rand_density: {:.2e}; nrand [1e6]: {:.2f}'.format(vol/1.e6, rand_density, nrand / 1.e6))

boundary_percent = 1.

zs      = np.arange(0.0, zmax+dz, dz)
Vs      = volcom(zs, Area) 

# https://docs.scipy.org/doc/scipy/reference/generated/scipy.interpolate.interp1d.html
Vz      = interp1d(Vs, zs, kind='linear', copy=True, bounds_error=True, fill_value=np.NaN, assume_sorted=False)

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

xyz      = cartesian(ras, decs, zs)

ras      = ras.astype(np.float32)
decs     = decs.astype(np.float32)
zs       = zs.astype(np.float32)
Vdraws   = Vdraws.astype(np.float32)
xyz      = xyz.astype(np.float32)

randoms           = Table(np.c_[ras, decs, zs, Vdraws], names=['RANDOM_RA', 'RANDOM_DEC', 'Z', 'V'])
randoms['RANDID'] = np.arange(len(randoms))
randoms['FIELD']  = gama_field(ras, decs)

assert  np.all(randoms['FIELD'].data == field)

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

randoms.meta = {'ZMIN': zmin,\
                'ZMAX': zmax,\
                'DZ': dz,\
                'NRAND': nrand,\
                'FIELD': field,\
                'Area': Area,\
                'BOUND_PERCENT': boundary_percent,\
                'VOL': vol,\
                'RAND_DENS': rand_density,\
                'VOL8': (4./3.)*np.pi*(8.**3.)}

randoms.meta['NRAND8']      = randoms.meta['VOL8'] * randoms.meta['RAND_DENS']
randoms.meta['NRAND8_PERR'] = np.sqrt(randoms.meta['NRAND8'])

print(randoms.meta)

runtime = calc_runtime(start, 'Writing {}'.format(opath), xx=randoms)

randoms.write(opath, format='fits', overwrite=True)

runtime = calc_runtime(start, 'Finished'.format(opath))
