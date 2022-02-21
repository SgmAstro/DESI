import os
import time
import numpy as np
import argparse
import astropy.io.fits   as     fits 

from   cosmo             import cosmo, volcom
from   scipy.interpolate import interp1d
from   astropy.table     import Table
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
realz   = args.realz

start   = time.time()

fields  = fetch_fields(survey)

assert  field in fields, f'Provided {field} field is not compatible with those available for {survey} survey ({fields})'

##  TODO: findfile.                                                                                                                                                                                  
##  opath = os.environ['RANDOMS_DIR'] + '/{}_{}_{:d}.fits'.format(prefix, field, realz)
opath   = findfile(ftype='randoms', dryrun=dryrun, field=field, survey=survey, prefix=prefix, realz=realz)

##  ras and decs.                                                                                                                                                              
if survey == 'gama':    
    from gama_limits import gama_field
    
    ra_min     = gama_limits[field]['ra_min']
    ra_max     = gama_limits[field]['ra_max']

    dec_min    = gama_limits[field]['dec_min']
    dec_max    = gama_limits[field]['dec_max']

    if dryrun == True:
        nrand  = 500
    
    else:
        nrand  = 1.e4

    ras        = ra_min
    ras[::2]   = ra_max
    decs       = np.arange(dec_min, dec_max, 1. / 60. / 60. / 10.) 
    
    np.random.shuffle(ras)
    np.random.shuffle(decs)

    randoms    = Table(np.c_[ras, decs], names=['RANDOM_RA', 'RANDOM_DEC'])

    ras        = np.arange(ra_min, ra_max, 1. / 60. / 60. / 10.)
    decs       = dec_min
    decs[::2]  = dec_max 

    np.random.shuffle(ras)
    np.random.shuffle(decs)

    randoms    = vstack((rand, Table(np.c_[ras, decs], names=['RANDOM_RA', 'RANDOM_DEC'])))
    nrand      = len(rand)

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


dz      = 1.e-6
zs      = np.arange(zmin, zmax+dz, dz)

np.random.shuffle(zs)

zs      = zs[:len(rand)]
Vs      = volcom(zs) - volcom(zmin)

print('Solved {:d} for field {}'.format(nrand, field))

randoms['Z']          = zs
randoms['V']          = volcom(zs) - volcom(zmin)
randoms['BOUNDID']    = np.arange(len(randoms))

randoms['FIELD']      = field
randoms['GAMA_FIELD'] = gama_field(ras, decs)

xyz                    = cartesian(ras, decs, zs)

randoms['CARTESIAN_X'] = xyz[:,0]
randoms['CARTESIAN_Y'] = xyz[:,1]
randoms['CARTESIAN_Z'] = xyz[:,2]

xyz                    = rotate(randoms['RANDOM_RA'], randoms['RANDOM_DEC'], xyz)

randoms['ROTCARTESIAN_X'] = xyz[:,0]
randoms['ROTCARTESIAN_Y'] = xyz[:,1]
randoms['ROTCARTESIAN_Z'] = xyz[:,2]

randoms.meta = {'ZMIN': zmin,\
                'ZMAX': zmax,\
                'DZ':     dz,\
                'NBOUND': nrand,\
                'FIELD': field,\
                'Area': Area}

print(randoms.meta)

runtime = calc_runtime(start, 'Writing {}'.format(opath), xx=randoms)

randoms.meta['EXTNAME'] = 'BOUNDARY'

randoms.write(opath, append=True)  

runtime = calc_runtime(start, 'Finished'.format(opath))
