import numpy as np
from cosmo import cosmo
from smith_kcorr import GAMA_KCorrection
from tmr_ecorr import tmr_ecorr, tmr_q

def theta(z, gmr):
    kcorr_r = GAMA_KCorrection(band='R')
    K = kcorr_r.k(z, gmr)
    E = tmr_ecorr(z, gmr, all=False)
    dist_mod = cosmo.distmod(z).value
    return dist_mod + K - E
