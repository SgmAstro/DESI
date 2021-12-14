import numpy as np

from   cosmo import cosmo
from   smith_kcorr import GAMA_KCorrection
from   tmr_ecorr import tmr_ecorr, tmr_q
from   cosmo import distmod

theta_kcorr_r = GAMA_KCorrection(band='R')
