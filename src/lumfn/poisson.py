import numpy as np


def fillfactor_expectation(bound_dist, radius=8.):
    # https://en.wikipedia.org/wiki/Spherical_cap/
    V_sphere    = (4./3.) * np.pi * radius ** 3.

    h           = np.abs(radius - bound_dist)
    h[bound_dist > radius] = 0.
    
    a           = np.sqrt(2 * h * radius - h * h)
    
    V_cap       = (np.pi * h / 6.) * (3. * a * a + h * h)
    V_res       = V_sphere - V_cap
    
    fill_factor = V_res / V_sphere
    
    return  fill_factor


def fillfactor_poisson(fill_factors, nbar=None, radius=None):
    if nbar == None:
        nbar = rand.meta['RAND_DENS']
    
    if radius == None:
        radius   = rand.meta['RSPHERE']
        V_sphere = rand.meta['VOL8']
        
    else:
        radius   = radius
        V_sphere = (4./3.) * np.pi * radius ** 3.
    
    NRAND8        = nbar * V_sphere * fill_factors
    SIGMA_NRAND8  = np.sqrt(NRAND8)
    SIGMA_FFACTOR = SIGMA_NRAND8 / NRAND8
    
    return  NRAND8, SIGMA_NRAND8, SIGMA_FFACTOR
