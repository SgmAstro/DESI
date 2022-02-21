import os
import time
import numpy as np
import argparse
import itertools
import astropy.io.fits   as     fits 

from   cosmo             import cosmo, volcom
from   scipy.interpolate import interp1d
from   astropy.table     import Table, vstack
from   cartesian         import cartesian, rotate
from   runtime           import calc_runtime
from   desi_randoms      import desi_randoms
from   findfile          import fetch_fields, findfile, overwrite_check
from   gama_limits       import gama_limits, gama_field


np.random.seed(314)

parser  = argparse.ArgumentParser(description='Calculate a set of boundary points')
parser.add_argument('-f', '--field',  type=str, help='select GAMA field [G9, G12, G15] or DESI rosette [R1...]', required=True)
parser.add_argument('-d', '--dryrun', help='Dryrun.', action='store_true')
parser.add_argument('-s', '--survey', help='Survey, e.g. GAMA, DESI, etc.', type=str, default='gama')
parser.add_argument('--sampling',     help='Sampling rate', default=90000)
parser.add_argument('--prefix',       help='filename prefix', default='randoms')
parser.add_argument('--nooverwrite',  help='Do not overwrite outputs if on disk', action='store_true')

# Defaults to GAMA Gold limits. 
parser.add_argument('--zmin', type=np.float32, help='Minimum redshift limit', default=0.039)
parser.add_argument('--zmax', type=np.float32, help='Maximum redshift limit', default=0.263)


args    = parser.parse_args()
field   = args.field.upper()
dryrun  = args.dryrun
survey  = args.survey.lower()
zmin    = args.zmin
zmax    = args.zmax
prefix  = args.prefix 
sampling = args.sampling
realz   = 0

start   = time.time()

fields  = fetch_fields(survey)

assert  field in fields, f'Provided {field} field is not compatible with those available for {survey} survey ({fields})'

##  TODO: findfile.                                                                                                                                                                                  
##  opath = os.environ['RANDOMS_DIR'] + '/{}_{}_{:d}.fits'.format(prefix, field, realz)
opath   = findfile(ftype='randoms', dryrun=dryrun, field=field, survey=survey, prefix=prefix, realz=realz)

##  ras and decs.                                                                                                                                                              
if survey == 'gama':    
    area       = 60. 

    ra_min     = gama_limits[field]['ra_min']
    ra_max     = gama_limits[field]['ra_max']

    dec_min    = gama_limits[field]['dec_min']
    dec_max    = gama_limits[field]['dec_max']

    pairs      = {'RA': (ra_min, ra_max), 'DEC': (dec_min, dec_max), 'Z': (zmin, zmax)}
    names      = list(pairs.keys())

    randoms    = []

    for key0 in names:
        keys         = list(pairs.keys())
        keys.remove(key0)
        
        key1         = keys[0]
        key2         = keys[1]

        print('Solving for {} boundary ({}, {})'.format(key0, key1, key2))
        
        pair0        = pairs[key0]
        pair1        = pairs[key1]
        pair2        = pairs[key2]

        continuous   = np.linspace(pair1[0], pair1[1], sampling)
        continuous   = np.tile(continuous, 2)

        np.random.shuffle(continuous)
         
        continuous2  = np.linspace(pair2[0], pair2[1], sampling)
        continuous2  = np.tile(continuous2, 2)

        np.random.shuffle(continuous2)

        discrete      = pair0[0] * np.ones_like(continuous)
        discrete[::2] = pair0[1]

        np.random.shuffle(discrete)

        to_add        = Table(np.c_[discrete, continuous, continuous2], names=['BOUND_{}'.format(key0), 'BOUND_{}'.format(key1), 'BOUND_{}'.format(key2)])
        to_add        = to_add['BOUND_RA', 'BOUND_DEC', 'BOUND_Z']

        randoms.append(to_add)

    randoms = vstack(randoms)
    randoms.rename_column('BOUND_Z', 'Z')

elif survey == 'desi':
    if 'NERSC_HOST' in os.environ.keys():
        raise NotImplementedError()

    else:
        print(f'As you are not running on nersc, the output of this script is assumed to be present at {opath} for dryrun: {dryrun}.')
        exit(0)

else:
    raise  NotImplementedError(f'No implementation for survey: {survey}')

if dryrun:
    nrand = 500
    opath = opath.replace('.fits', '_dryrun.fits')

else:
    nrand = len(randoms)

print('Solved {:d} for field {}'.format(nrand, field))

randoms.pprint()

randoms['V']          = volcom(randoms['Z'], area=area) - volcom(zmin, area=area)
randoms['BOUNDID']    = np.arange(len(randoms))

randoms['FIELD']      = field
randoms['GAMA_FIELD'] = gama_field(randoms['BOUND_RA'], randoms['BOUND_DEC'])

xyz                    = cartesian(randoms['BOUND_RA'], randoms['BOUND_DEC'], randoms['Z'])

randoms['CARTESIAN_X'] = xyz[:,0]
randoms['CARTESIAN_Y'] = xyz[:,1]
randoms['CARTESIAN_Z'] = xyz[:,2]

xyz                    = rotate(randoms['BOUND_RA'], randoms['BOUND_DEC'], xyz)

randoms['ROTCARTESIAN_X'] = xyz[:,0]
randoms['ROTCARTESIAN_Y'] = xyz[:,1]
randoms['ROTCARTESIAN_Z'] = xyz[:,2]

randoms.meta = {'ZMIN': zmin,\
                'ZMAX': zmax,\
                'NBOUND': nrand,\
                'FIELD': field,\
                'SAMPLING': sampling,\
                'Area': area}

print(randoms.meta)

randoms.meta['EXTNAME'] = 'BOUNDARY'

if os.path.isfile(opath):
    runtime = calc_runtime(start, f'Appending BOUNDARY extension to {opath}', xx=randoms)

else:
    raise  RuntimeError(f'Failed to find {opath} needed to append.')

randoms.write(opath, append=True, overwrite=True)  

runtime = calc_runtime(start, 'Finished'.format(opath))
