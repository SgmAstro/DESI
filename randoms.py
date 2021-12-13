import numpy as np

from cosmo import cosmo
from scipy.interpolate import interp1d


np.random.seed(314)

dz   = 1.e-4

zmin = 300. / 2.9979e5
zmax = 0.6

zs   = np.arange(zmin, zmax, dz)
Vs   = cosmo.comoving_volume(zs) * u.Mpc**-3

# https://docs.scipy.org/doc/scipy/reference/generated/scipy.interpolate.interp1d.html
Vz   = interp1d(zs, Vs, kind='linear', copy=True, bounds_error=True, fill_value=nan, assume_sorted=True)
