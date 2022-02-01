import os
import numpy as np

from   findfile import findfile


def lum_binner(x, dM=0.1):
    '''
    Eqn. 2.10a, W(x), of Efstathiou, Ellis & Peterson.
    '''
    return  np.abs(x) <= (dM / 2.)

def lum_visible(x, dM=0.1):
    '''                                                                                                                          
    Eqn. 2.10b, H(x), of Efstathiou, Ellis & Peterson.                                                                            
    '''
    if x >= (dM/2.):
        return  0.0

    elif x <= (dM/2.):
        return  1.0

    else:
        return  -x/dM + 1./2.

def lumfn_stepwise_eval(vmax, phis, phi_Ms, Mcol='MCOLOR_0P0', Mmin_col='DDPMALL_0P0_VISZ'):
    '''
    Eqn. 2.12, of Efstathiou, Ellis & Peterson.   
    '''
    Ms    = vmax[Mcol]
    Mmins = vmax[Mmin_col]

    dM    = np.diff(phi_Ms)
    dM    = dM[0]
    
    num   = np.count_nonzero(lum_binner(Ms - phi_Ms))

    facs  = []

    for Mmin in Mmins:
        fac   = lum_visible(phi_Ms - Mmin)
        fac  *= dM * phis
        fac   = np.sum(fac)
        
        facs.append(fac)

    fac   = np.array(fac)
        
    den   = lum_visible(Ms - Mmins)
    den  /= fac
    
    den   = np.sum(den)

    # dM * phis.
    return  num / den
    
def lumfn_stepwise(vmax, Mcol='MCOLOR_0P0', Mmin_col='DDPMALL_0P0_VISZ'):
    phi_Ms   = np.arange(-26., -16., 0.5) 
    phi_init = np.ones_like(phi_Ms)

    diff     = 1.e99
    
    while  (diff > 1.e-6):
        new_phis = lumfn_stepwise_eval(vmax, phi_init, phi_Ms, Mcol=Mcol, Mmin_col=Mmin_col)
        diff     = (new_phis - phi_init)**2.

        #  Update previous estimate. 
        phi_init = new_phis
        
    return  phi_Ms, new_phis
        

if __name__ == '__main__':
    from astropy.table import Table


    print('\n\n')

    fpath = findfile('vmax')

    print(fpath)

    vmax = Table.read(findfile('vmax'))
    vmax.pprint()

    result = lumfn_stepwise(vmax)

    print('\n\nDone.\n\n')
