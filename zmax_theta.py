import numpy as np

from   cosmo import cosmo
from   smith_kcorr import GAMA_KCorrection
from   tmr_ecorr import tmr_ecorr, tmr_q


def theta(zs, rest_gmrs, all=False):
    theta_kcorr_r = GAMA_KCorrection(band='R')

    theta = []
    
    for z, rest_gmr in zip(zs, rest_gmrs):
        K = theta_kcorr_r.k(z, rest_gmr)
        E = tmr_ecorr(z, rest_gmr, all=all)
        
        dist_mod = cosmo.distmod(z)
        
        thetas.append(dist_mod + K - E)

    return thetas
