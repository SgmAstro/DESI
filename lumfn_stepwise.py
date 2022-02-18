import os
import numpy as np

from   findfile import findfile
from   schechter import named_schechter


def lum_binner(x, dM=0.1):
    '''
    Eqn. 2.10a, W(x), of Efstathiou, Ellis & Peterson.
    '''
    return  np.abs(x) <= (dM / 2.)

def lum_visible(x, dM=0.1):
    '''                                                                                                                          
    Eqn. 2.10b, H(x), of Efstathiou, Ellis & Peterson.                                                                            
    '''

    result = -x/dM + 1./2.

    result[x >=  (dM / 2.)] = 0.0
    result[x <= -(dM / 2.)] = 1.0
    
    return  result

def lumfn_stepwise_eval(vmax, phi_M, phi, phis, phi_Ms, dM, Mcol='DDPMALL_0P0', Mmin_col='DDPMALL_0P0_VISZ', vmax_stepwise=False):
    '''
    Eqn. 2.12, of Efstathiou, Ellis & Peterson.   
    '''
    Ms    = vmax[Mcol]
    Mmins = vmax[Mmin_col]

    num   = np.count_nonzero(lum_binner(phi_M - Ms))

    facs  = lum_visible(phi_M - Mmins) 
    facs[facs > 0.0] /= np.sum(dM * phis[:,None] * lum_visible(phi_Ms[:,None] - Mmins[facs > 0.0]), axis=0)

    den   = np.sum(facs) 

    if vmax_stepwise:
        return  den
     
    # dM * phis.
    return  num / den
    
def lumfn_stepwise(vmax, Mcol='DDPMALL_0P0', Mmin_col='DDPMALL_0P0_VISZ'):
    dM        = 0.1

    phi_Ms    = np.arange(-26., -16., dM) 
    phi_init  = 1.e-2 * np.ones_like(phi_Ms)

    diff      = 1.e99
    phis      = phi_init 

    iteration = 0

    while  (diff > 1.e-6):
        print('Solving for iteration {:d} with diff. {:.6e}'.format(iteration, diff))
        
        new_phis = []
    
        for i, (phi_M, phi) in enumerate(zip(phi_Ms, phis)):
            new_phis.append(lumfn_stepwise_eval(vmax, phi_M, phi, phis, phi_Ms, dM, Mcol=Mcol, Mmin_col=Mmin_col))

        new_phis   = np.array(new_phis)
        
        diff       = np.sum((new_phis - phis)**2.)
     
        #  Update previous estimate. 
        phis       = new_phis

        iteration += 1

    _vmaxs = []

    for i, (phi_M, phi) in enumerate(zip(phi_Ms, phis)):
        _vmaxs.append(lumfn_stepwise_eval(vmax, phi_M, phi, phis, phi_Ms, dM, Mcol=Mcol, Mmin_col=Mmin_col, vmax_stepwise=True))

    vmaxs_stepwise = np.array(_vmaxs)

    phis  = phis / dM 

    nn    = dM * np.sum(phis)

    phis /= nn

    return  phi_Ms, phis, vmaxs_stepwise
        

if __name__ == '__main__':
    import pylab as pl
    
    from   astropy.table import Table


    print('\n\n')

    fpath = findfile('vmax')

    print(fpath)

    fpath = '/cosma5/data/durham/dc-moor2/GAMA4/gama_gold_ddp.fits'

    vmax = Table.read(fpath)

    vmax.pprint()

    phi_Ms, phis, vmaxs_stepwise = lumfn_stepwise(vmax)

    np.savetxt('stepwise.txt', np.c_[phi_Ms, phis, vmaxs_stepwise])
    
    print('\n\nDone.\n\n')
