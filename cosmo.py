import numpy as np
import astropy.units as u

from    astropy.cosmology       import  FlatLambdaCDM

# setting cosmological parameters
h = 1
cosmo = FlatLambdaCDM(H0=100*h * u.km / u.s / u.Mpc, Tcmb0=2.725 * u.K, Om0= 0.25)

def distmod(zs):
    return 5. * np.log10(cosmo.luminosity_distance(zs).value) - 25.
