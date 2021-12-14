import os
import time
import numpy as np

from cosmo import distmod
from smith_kcorr import GAMA_KCorrection
from tmr_ecorr import tmr_ecorr
from scipy.optimize import brentq, minimize
from astropy.table import Table

kcorr_r = GAMA_KCorrection(band='R')

def theta(z, rest_gmr_0p1, rest_gmr_0p0, aall=True):
    z = np.atleast_1d(z)
    rest_gmr_0p1 = np.atleast_1d(rest_gmr_0p1)
    rest_gmr_0p0 = np.atleast_1d(rest_gmr_0p0)
    
    result = distmod(z) + kcorr_r.k_nonnative_zref(0.0, z, rest_gmr_0p1) + tmr_ecorr(z, rest_gmr_0p0, aall=aall)

    return result[0]
    
def solve_theta(rest_gmr_0p1, rest_gmr_0p0, thetaz, dr, aall=True):
     def diff(x):
          return theta(x, rest_gmr_0p1, rest_gmr_0p0, aall=aall) - thetaz - dr

     def absdiff(x):
        return np.abs(diff(x))

     warn = 0

     try:
        result = brentq(diff, 1.e-3, 0.6)

     except ValueError as VE:
        warn = 1

        # Brent method fails, requires sign change across boundaries.                                                                                          
        result = minimize(absdiff, 0.3)

        if result.success:
            result = result.x[0]

        else:
             warn = 2
             result = -99.

     return  result, warn

def zmax(rest_gmrs_0p1, rest_gmrs_0p0, theta_zs, drs, aall=True, debug=True):
   result = []
   start = time.time()

   if debug:
        print('Solving for zlimit.')

   for i, (rest_gmr_0p1, rest_gmr_0p0, theta_z, dr) in enumerate(zip(rest_gmrs_0p1, rest_gmrs_0p0, theta_zs, drs)):
        interim, warn = solve_theta(rest_gmr_0p1, rest_gmr_0p0, theta_z, dr, aall=aall)

        result.append([interim, warn])

        if (i % 500 == 0) & debug:
             runtime = (time.time() - start) / 60.

             print('{:.3f}% complete after {:.2f} mins.'.format(100. * i / len(theta_zs), runtime))

   result = np.array(result)

   return  result[:,0], result[:,1]


ngal=5000
#ngal=20000
nproc=4

rlim = 19.8
rmax = 12.0

root = os.environ['CSCRATCH'] + '/norberg/'
fpath = root + '/GAMA4/gama_gold_kE_{:d}k.fits'.format(np.int(ngal / 1000.))

dat = Table.read(fpath)
dat.pprint()

dat['DELTA_RPETRO_FAINT'] = rlim - dat['R_PETRO']

zmaxs, warn = zmax(dat['REST_GMR_0P1'], dat['REST_GMR_0P0'], dat['Z_THETA_QCOLOR'], dat['DELTA_RPETRO_FAINT'], aall=True, debug=True)

dat['ZMAX'] = zmaxs
dat['ZMAX_WARN'] = warn

dat['DELTA_RPETRO_BRIGHT'] = rmax - dat['R_PETRO']

zmins, warn = zmax(dat['REST_GMR_0P1'], dat['REST_GMR_0P0'], dat['Z_THETA_QCOLOR'], dat['DELTA_RPETRO_BRIGHT'], aall=True, debug=True)

dat['ZMIN'] = zmins
dat['ZMIN_WARN'] = warn

dat.write(root + '/GAMA4/gama_gold_zmax_{:d}k.fits'.format(np.int(ngal / 1000.)), format='fits', overwrite=True)
