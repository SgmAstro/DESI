import os
import sys
import time
import numpy as np
import argparse

from   cosmo             import cosmo, volcom
from   scipy.interpolate import interp1d
from   astropy.table     import Table
from   cartesian         import cartesian, rotate
from   runtime           import calc_runtime
from   desi_randoms      import desi_randoms
from   findfile          import fetch_fields, findfile, overwrite_check, call_signature
from   gama_limits       import gama_limits, gama_field
from   bitmask           import lumfn_mask, consv_mask

def randoms(field='G9', survey='gama', density=1., zmin=0.039, zmax=0.263, dryrun=False, prefix='', seed=314, oversample=8, realz=0):
    start   = time.time()

    fields  = fetch_fields(survey)

    assert  field in fields, f'Provided {field} field is not compatible with those available for {survey} survey ({fields})'

    opath   = findfile(ftype='randoms', dryrun=dryrun, field=field, survey=survey, prefix=prefix, realz=realz, oversample=oversample)

    if args.nooverwrite:
        overwrite_check(opath)

    seed    = seed + realz + 50 * oversample

    np.random.seed(seed)

    call_signature(dryrun, sys.argv)

    ##  ras and decs.                                                                                                                                                              
    if survey == 'gama':    
        from gama_limits import gama_field

        Area    = 60.

        # TODO:
        # field   = fields[args.field]

        ra_min  = gama_limits[field]['ra_min']
        ra_max  = gama_limits[field]['ra_max']

        dec_min = gama_limits[field]['dec_min']
        dec_max = gama_limits[field]['dec_max']

        ctheta_min = np.cos(np.pi/2. - np.radians(dec_min))
        ctheta_max = np.cos(np.pi/2  - np.radians(dec_max))

        vol          = volcom(zmax, Area) - volcom(zmin, Area)

        if dryrun == True:
            ndryrun = 2000
            nrand   = ndryrun

        else:
            nrand     = int(np.ceil(vol * density * oversample) / 2.0)
            print('NRAND IS:', nrand)

        cos_theta = np.random.uniform(ctheta_min, ctheta_max, nrand)
        theta     = np.arccos(cos_theta)
        decs      = np.pi/2. - theta
        decs      = np.degrees(decs)

        ras       = np.random.uniform(ra_min, ra_max, nrand)

        randoms   = Table(np.c_[ras, decs], names=['RANDOM_RA', 'RANDOM_DEC'])

    elif survey == 'desi':
        if 'NERSC_HOST' in os.environ.keys():
            # Support to run on nersc only.
            randoms = desi_randoms(ros=int(field[1:]))
            nrand   = len(randoms)

            # TODO: add dryrun nrand fix (as above in GAMA)

            # Original density of 2500 per sq. deg. 
            Area    = nrand / 2500. 

        elif 'ddp1' in prefix:
            # Assume you are on cosma, rewriting ddp1-like redshift limits to ddp1 randoms based on the assumed present randoms. 
            randoms = findfile(ftype='randoms', dryrun=dryrun, field=field, survey=survey, prefix=None, realz=realz)

            print(f'As you are not running on nersc, an input of this script is assumed to be present at {randoms} for dryrun: {dryrun}.')

            randoms = Table.read(randoms)
            nrand   = len(randoms)

            Area    = nrand / 2500.

        else:
            print(f'As you are not running on nersc, the output of this script is assumed to be present at {opath} for dryrun: {dryrun}.')
            exit(0)

    else:
        raise  NotImplementedError(f'No implementation for survey: {survey}')
    
    ##  Vs and zs.
    dz      = 1.e-4

    Vmin    = volcom(zmin, Area)
    Vmax    = volcom(zmax, Area)

    vol     = Vmax - Vmin

    density = nrand / vol

    if dryrun:
        # TODO: Hard coded above, don't hard code twice. 
        nrand = ndryrun

    if not os.path.isdir(os.environ['RANDOMS_DIR']):
        print('Creating {}'.format(os.environ['RANDOMS_DIR']))

        os.makedirs(os.environ['RANDOMS_DIR'])

    print('Volume [1e6]: {:.2f}; oversample: {:.2f};  density: {:.2e}; nrand [1e6]: {:.2f}'.format(vol/1.e6, oversample, density, nrand / 1.e6))

    zs      = np.arange(0.0, zmax+dz, dz)
    Vs      = volcom(zs, Area) 

    # https://docs.scipy.org/doc/scipy/reference/generated/scipy.interpolate.interp1d.html
    Vz      = interp1d(Vs, zs, kind='linear', copy=True, bounds_error=True, fill_value=np.NaN, assume_sorted=False)

    Vdraws  = np.random.uniform(0., 1., nrand)
    Vdraws  = Vmin + Vdraws * (Vmax - Vmin)

    zs      = Vz(Vdraws)

    print('Solved {:d} for field {}'.format(nrand, field))

    print('Applying rotation.')

    ras      = randoms['RANDOM_RA']
    decs     = randoms['RANDOM_DEC']

    print(len(ras), len(decs), len(zs))

    xyz      = cartesian(ras, decs, zs)
    
    '''
    ras      = ras.astype(np.float32)
    decs     = decs.astype(np.float32)
    zs       = zs.astype(np.float32)
    Vdraws   = Vdraws.astype(np.float32)
    xyz      = xyz.astype(np.float32)
    '''

    randoms['Z'] = zs
    randoms['V'] = Vdraws
    randoms['RANDID'] = np.arange(len(randoms))

    randoms['FIELD']      = field
    randoms['GAMA_FIELD'] = gama_field(ras, decs)

    # assert  np.all(randoms['FIELD'].data == field)

    randoms['CARTESIAN_X'] = xyz[:,0]
    randoms['CARTESIAN_Y'] = xyz[:,1]
    randoms['CARTESIAN_Z'] = xyz[:,2]

    xyz = rotate(randoms['RANDOM_RA'], randoms['RANDOM_DEC'], xyz)

    randoms['ROTCARTESIAN_X'] = xyz[:,0]
    randoms['ROTCARTESIAN_Y'] = xyz[:,1]
    randoms['ROTCARTESIAN_Z'] = xyz[:,2]

    '''
    elif survey == 'desi':    
        randoms['IS_BOUNDARY'][randoms['ROS_DIST']   > np.percentile(randoms['ROS_DIST'],   100. - boundary_percent)] = 1
        randoms['IS_BOUNDARY'][randoms['ROS_DIST']   < np.percentile(randoms['ROS_DIST'],   boundary_percent)]        = 1
    '''

    randoms['ZSURV']        = randoms['Z']
    randoms['IN_D8LUMFN']   = np.zeros_like(randoms['FIELD'], dtype=int)
    randoms['CONSERVATIVE'] = np.zeros_like(randoms['FIELD'], dtype=int)

    randoms.meta = {'ZMIN':   zmin,\
                    'ZMAX':   zmax,\
                    'DZ':       dz,\
                    'NRAND': nrand,\
                    'FIELD': field,\
                    'Area':   Area,\
                    'VOL':     vol,\
                    'RAND_DENS': density,\
                    'VOL8': (4./3.)*np.pi*(8.**3.)}

    randoms.meta['NRAND8']      = randoms.meta['VOL8'] * randoms.meta['RAND_DENS']
    randoms.meta['NRAND8_PERR'] = np.sqrt(randoms.meta['NRAND8'])

    print(randoms.meta)

    runtime = calc_runtime(start, 'Writing {}'.format(opath), xx=randoms)

    randoms.write(opath, format='fits', overwrite=True)

    runtime = calc_runtime(start, 'Finished'.format(opath))



if __name__ == '__main__':
    parser  = argparse.ArgumentParser(description='Select GAMA field.')
    parser.add_argument('-f', '--field',  type=str, help='select GAMA field [G9, G12, G15] or DESI rosette [R1...]', default='G9')
    parser.add_argument('-d', '--dryrun', help='Dryrun.', action='store_true')
    parser.add_argument('-s', '--survey', help='Survey, e.g. GAMA, DESI, etc.', type=str, default='gama')
    parser.add_argument('--realz',        help='Realization', default=0, type=int)
    parser.add_argument('--prefix',       help='filename prefix', default='randoms')
    parser.add_argument('--nooverwrite',  help='Do not overwrite outputs if on disk', action='store_true')
    parser.add_argument('--density',      help='Random density per (Mpc/h)^3', default=1.0, type=float)
    parser.add_argument('--oversample',   help='Oversampling factor for fillfactor counting.', default=8, type=int)
    parser.add_argument('--seed',         help='Random seed.', default=314, type=int)
    
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
    seed    = args.seed

    density = args.density
    oversample = args.oversample

    assert oversample in np.arange(1, 21, 1)
    
    for xx in [1, oversample]:        
        randoms(field=field, survey=survey, density=density, zmin=zmin, zmax=zmax, dryrun=dryrun, prefix=prefix, seed=seed, oversample=xx, realz=realz)


