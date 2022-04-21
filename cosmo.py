import numpy as np
import astropy.units as u

from   astropy.cosmology       import  FlatLambdaCDM

# setting cosmological parameters
h     = 1
cosmo = FlatLambdaCDM(H0=100*h * u.km / u.s / u.Mpc, Tcmb0=2.725 * u.K, Om0= 0.25)

def fsky(area_sqdeg):
    return area_sqdeg / (360.**2. / np.pi)

def negz_proof(func):
    '''
    Decorator to handle negative redshifts.
    '''
  
    def wrap(zs):
        negz     = zs <= 0.0
        
        zs[negz] = 1.e-99

        result   = func(zs)
        result[negz] = np.nan
 
        return result

    return wrap

@negz_proof
def distmod(zs):
    return 5. * np.log10(cosmo.luminosity_distance(zs).value) + 25.

@negz_proof
def distcom(zs):
    return cosmo.comoving_distance(zs).value

@negz_proof
def volcom(zs, area):
    return (4./3.) * np.pi * fsky(area) * distcom(zs)**3.


if __name__ == '__main__':
    zs  = np.arange(-10., 10., 1.)
    mus = distmod(zs)

    for z, mu in zip(zs, mus):
        print('{:.6f}\t{:.6f}'.format(z, mu))
