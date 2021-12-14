import numpy as np

from   cosmo import cosmo
from   smith_kcorr import GAMA_KCorrection
from   tmr_ecorr import tmr_ecorr, tmr_q
from   cosmo import distmod


def zmax_theta(zs, rest_gmrs_0p1, rest_gmrs_0p0, aall=False):
    theta_kcorr_r = GAMA_KCorrection(band='R')
        
    thetas = []

    for z, rest_gmr_0p1, rest_gmr_0p0 in zip(zs, rest_gmrs_0p1, rest_gmrs_0p0):
        z = np.atleast_1d(z)
        rest_gmr_0p1 = np.atleast_1d(rest_gmr_0p1)
        rest_gmr_0p0 = np.atleast_1d(rest_gmr_0p0)
        
        K = theta_kcorr_r.k(z, rest_gmr_0p1)
        E = tmr_ecorr(z, rest_gmr_0p0, aall=aall)
        
        dist_mod = distmod(z)
        
        thetas.append(dist_mod + K + E)

    thetas = np.array(thetas)
        
    return thetas
