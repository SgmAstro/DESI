import os
import sys
import tqdm
import time
import argparse
import runtime
import numpy           as     np
import multiprocessing

from   cosmo           import distmod, volcom
from   smith_kcorr     import GAMA_KCorrection
from   tmr_ecorr       import tmr_ecorr
from   scipy.optimize  import brent, minimize, brentq
from   astropy.table   import Table
from   functools       import partial
from   multiprocessing import Pool
from   findfile        import findfile, overwrite_check, write_desitable, fetch_header
from   config          import Configuration
from   abs_mag         import abs_mag


kcorr_r          = GAMA_KCorrection(band='R')

def theta(z, rest_gmr_0p1, rest_gmr_0p0, thetaz=None, dr=None, aall=False, absolute=False):
    z            = np.atleast_1d(z)
    rest_gmr_0p1 = np.atleast_1d(rest_gmr_0p1)
    rest_gmr_0p0 = np.atleast_1d(rest_gmr_0p0)
    
    result       = distmod(z) + kcorr_r.k_nonnative_zref(0.0, z, rest_gmr_0p1) + tmr_ecorr(z, rest_gmr_0p0, aall=aall)

    if thetaz != None:
        result -= thetaz
        result -= dr

    if absolute:
        result  = np.abs(result) 

    return  result[0]
    
def solve_theta(rest_gmr_0p1, rest_gmr_0p0, thetaz, dr, aall=False, debug=False, startz=None):
    if startz == None:
        startz = 2.5

    try:
        result = brentq(theta, 1.e-6, 1.6, args=(rest_gmr_0p1, rest_gmr_0p0, thetaz, dr, aall, False))
        warn   = 0
        method = 0

    except ValueError as VE:
        if debug:
            print(VE)

        # Brent method fails, requires sign change across boundaries.                                                                                          
        result = minimize(theta, startz, args=(rest_gmr_0p1, rest_gmr_0p0, thetaz, dr, aall, True), method='Nelder-Mead')

        if result.success:
            result = result.x[0]
            warn   = 0
            method = 1

        else:
             try:
                 result = brent(theta, brack=(1.e-6, 1.6), args=(rest_gmr_0p1, rest_gmr_0p0, thetaz, dr, aall, True))
                 warn   = 0
                 method = 2

             except ValueError as VE:
                 result = -99
                 warn   =   1 

    return  result, warn, method

def zmax(rest_gmrs_0p1, rest_gmrs_0p0, theta_zs, drs, aall=False, debug=True, nproc=14, startz=None):
   result = []
   start  = time.time()

   if debug:
        print('Solving for zlimit.')

   with multiprocessing.get_context('spawn').Pool(processes=nproc) as pool:
       arglist = list(zip(rest_gmrs_0p1, rest_gmrs_0p0, theta_zs, drs))
       result  = pool.starmap(partial(solve_theta, aall=aall, startz=startz), arglist)
   
       pool.close()

       # https://stackoverflow.com/questions/38271547/when-should-we-call-multiprocessing-pool-join                                                                                                        
       pool.join()

   result = np.array(result)

   return  result[:,0], result[:,1], result[:,2]


if __name__ == '__main__':
    parser  = argparse.ArgumentParser(description='Gen zmax cat.')
    parser.add_argument('--log', help='Create a log file of stdout.', action='store_true')  
    parser.add_argument('-a', '--aall',   help='All Q, no red/blue split.', action='store_true')
    parser.add_argument('-d', '--dryrun', help='Dryrun.', action='store_true')
    parser.add_argument('--nproc', type=int, help='Number of processors', default=14)
    parser.add_argument('-s', '--survey', help='Select survey', default='gama')
    parser.add_argument('--theta_def',    help='Specifier for definition of theta', default='Z_THETA_QCOLOR')
    parser.add_argument('--config',       help='Path to configuration file', type=str, default=findfile('config'))
    parser.add_argument('--nooverwrite',  help='Do not overwrite outputs if on disk', action='store_true')
    
    args      = parser.parse_args()
    log       = args.log
    aall      = args.aall
    nproc     = args.nproc 
    dryrun    = args.dryrun
    survey    = args.survey.lower()
    theta_def = args.theta_def

    config    = Configuration(args.config)
    config.update_attributes('zmax', args)
    config.write()

    rlim      = fetch_header(ftype='gold', name='RLIM', survey=survey)
    rmax      = fetch_header(ftype='gold', name='RMAX', survey=survey)

    start     = time.time()

    if log:
        logfile = findfile(ftype='zmax', dryrun=False, survey=survey, log=True)

        print(f'Logging to {logfile}')

        sys.stdout = open(logfile, 'w')

    print('Assuming {:.4f} < r < {:.4f}'.format(rmax, rlim))
    print('Assuming Q ALL = {}'.format(aall))
    
    fpath  = findfile(ftype='kE',   dryrun=dryrun, survey=survey)
    opath  = findfile(ftype='zmax', dryrun=dryrun, survey=survey)
    
    if args.nooverwrite:
        overwrite_check(opath)

    print('Reading {}.'.format(fpath))
        
    dat = Table.read(fpath)
    dat.pprint()

    dat['DELTA_DETMAG_FAINT'] = rlim - dat['DETMAG']

    print('Solving for {} bounding curve'.format(rlim))
    
    zmaxs, warn, method = zmax(dat['REST_GMR_0P1'],\
                               dat['REST_GMR_0P0'],\
                               dat[theta_def],\
                               dat['DELTA_DETMAG_FAINT'],\
                               aall=aall,\
                               nproc=nproc,\
                               debug=True)

    dat['ZMAX']           = zmaxs
    dat['ZMAX_WARN']      = warn
    dat['ZMAX_METHOD']    = method

    print('Solving for {} bounding curve'.format(rmax))

    dat['DELTA_DETMAG_BRIGHT'] = rmax - dat['DETMAG']
    
    zmins, warn, method = zmax(dat['REST_GMR_0P1'],\
                               dat['REST_GMR_0P0'],\
                               dat[theta_def],\
                               dat['DELTA_DETMAG_BRIGHT'],\
                               aall=aall,\
                               nproc=nproc,\
                               startz=0.1,\
                               debug=True)

    dat['ZMIN']           = zmins
    dat['ZMIN_WARN']      = warn
    dat['ZMIN_METHOD']    = method

    dat.meta['THETA_DEF'] = theta_def

    dat['VMAX']  = volcom(dat['ZMAX'], dat.meta['AREA'])
    dat['VMAX'] -= volcom(dat['ZMIN'], dat.meta['AREA'])

    print('Writing {}.'.format(opath))

    dat.pprint()

    write_desitable(opath, dat)

    nwarn   = (dat['ZMAX_WARN'].data > 0) | (dat['ZMIN_WARN'].data > 0)
    nwarn   = np.count_nonzero(nwarn)

    towarn  = dat[(dat['ZMAX_WARN'].data > 0) | (dat['ZMIN_WARN'].data > 0)]

    if len(towarn) > 0:
        print(f'WARNING:  zmax/min warnings triggered on {nwarn} galaxies.')

        towarn.pprint()


    print(f'\n\nMETHOD\tZMAX\tZMIN')

    for method, name in zip(range(3), ['BRENTQ', 'NELDER', 'BRENT']):
        print('{}\t{:d}\t{:d}'.format(name, np.count_nonzero(dat['ZMAX_METHOD'] == method), np.count_nonzero(dat['ZMIN_METHOD'] == method)))


    runtime = (time.time() - start) / 60.

    print('\n\nDone in {} mins.\n\n'.format(runtime))

    if log:
        sys.stdout.close()
