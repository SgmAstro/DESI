import numpy as np
import scipy.integrate as integrate

from   data.schechters import schechters


def schechter(M, phistar, Mstar, alpha):
    expa         = 10. ** (0.4 * (Mstar - M) * (1. + alpha))
    expb         = np.exp(-10. ** (0.4 * (Mstar - M)))

    return  np.log(10.) * phistar * expa * expb / 2.5

def schechter_d8(M, d8, params=False, fit=True):
    '''
    TMR d8 Schechter.
    '''

    alpha         = -1.25
    Mstar         = -20.70 - 0.67 * np.log10(1. + d8) 
    log10phistar  = -2.030 + 1.01 * np.log10(1. + d8)

    if not fit:
        # Middle panel TMR Fig. 7 shows least dense bin is above fit
        # by approx. 0.3 dex.
        if d8 < -0.75:
            alpha         = -1.375
            log10phistar +=  0.300

    phistar       = 10.**log10phistar

    if params:
        return  np.log10(1. + d8), log10phistar, Mstar, alpha  

    else:
        return  schechter(M, phistar, Mstar, alpha)

def named_schechter(M, named_type='TMR', zz=None, evolve=False):
    params       = schechters[named_type]
    
    log10phistar = params['log10phistar']
    Mstar        = params['Mstar']
    alpha        = params['alpha'] 

    P            = params['P']
    Q            = params['Q']
    zref         = params['zref']

    if zz == None:
        zz       = zref

    zz           = np.array([zz], copy=True)[0]
    
    phistar      = 10. ** log10phistar

    if evolve:
        Mstar       -= Q * (zz - zref)
        phistar     *= 10. ** (0.4 * P * (zz - zref))

    return schechter(M, phistar, Mstar, alpha)
