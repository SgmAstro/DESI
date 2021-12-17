import numpy             as np
import matplotlib.pyplot as plt
import cosmo             as cosmo
import os

from   astropy.table     import Table
from   smith_kcorr       import GAMA_KCorrection, GAMA_KCorrection_color
from   rest_gmr          import smith_rest_gmr
from   tmr_ecorr         import tmr_ecorr, tmr_q
from   abs_mag           import abs_mag
from   data.ke_params    import *


kcorr_r  = GAMA_KCorrection(band='R')
kcorr_RG = GAMA_KCorrection_color()

# To be looped over for a total of 7 x 3 x 2 curves.
gmrs_0p1 =  np.array([0.131, 0.298, 0.443, 0.603, 0.785, 0.933, 1.067])  
gmrs_0p0 =  np.array([0.158, 0.298, 0.419, 0.553, 0.708, 0.796, 0.960])

rlims    = [12., 19.8] # bright and faint limits.

zs       = np.arange(0.01, 0.6, 0.01)
mus      = cosmo.distmod(zs)

root     = os.environ['CSCRATCH'] + '/norberg/GAMA4/ddrp_limits/'

count    = 0

for rlim in rlims:
    rs = rlim * np.ones_like(zs)

    for aall, all_type in zip([True, False], ['QALL', 'QCOLOR']):
        for gmr_0P1 in gmrs_0p1:
            gmr_0P1  = gmr_0P1 * np.ones_like(zs)
            gmr_0P0  = kcorr_RG.rest_gmr_nonnative(gmr_0P1)

            ks       = kcorr_r.k_nonnative_zref(0.0, zs, gmr_0P1)
            es       = tmr_ecorr(zs, gmr_0P0, aall=aall)
            Mrs_0P0  = abs_mag(rs, mus, ks, es)

            dat      = Table(np.c_[zs, ks, es, Mrs_0P0], names=['Z', 'K', 'E', 'M0P0_{}'.format(all_type)])
            
            opath    = root + 'ddrp_limit_{:d}.fits'.format(count)
            
            dat.meta = {'RLIM': rlim, 'ALL': aall, 'GMR_0P1': gmr_0P1[0], 'GMR_0P0': gmr_0P0[0]}
            
            dat.write(opath, format='fits', overwrite=True)
            
            count   += 1

            print('Solved for {} {} {}'.format(rlim, all_type, gmr_0P1))
