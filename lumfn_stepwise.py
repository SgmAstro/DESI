import  os
import  tqdm
import  argparse
import  numpy           as     np

from    multiprocessing import Pool
from    astropy.table   import Table
from    findfile        import findfile, overwrite_check
from    schechter       import named_schechter
from    ddp             import initialise_ddplimits


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

def lumfn_stepwise_eval(vmax, phi_M, phi, phis, phi_Ms, dM, Mcol='MALL_0P0', Mmin_col='DDPMALL_0P0_VISZ', survey='gama'):
    '''
    Eqn. 2.12, of Efstathiou, Ellis & Peterson.   
    '''
    Ms    = vmax[Mcol]
    Mmins = vmax[Mmin_col]

    num   = np.count_nonzero(lum_binner(phi_M - Ms))

    bright_curve, bright_curve_r, faint_curve, faint_curve_r = initialise_ddplimits(survey=survey)

    # HACK MALL, MQALL?
    assert  Mcol == 'MALL_0P0' 

    zmin  = bright_curve(phi_M) 
    zmax  = faint_curve(phi_M)

    '''
    with Pool(nproc) as pool:
        result = pool.imap(process_one, splits)                                                                                                                                                            
    '''

    # Deprecated:
    # facs              = lum_visible(phi_M - Mmins) 
    # facs[facs > 0.0] /= np.sum(dM * phis[:,None] * lum_visible(phi_Ms[:,None] - Mmins[facs > 0.0]), axis=0)
    # den               = np.sum(facs) 

    Mmin    = bright_curve(vmax['Z'])
    Mmax    =  faint_curve(vmax['Z'])
    
    vol_lim = vmax[(vmax['Z'] > zmin) * (vmax['Z'] < zmax)]

    # Ngal x NM;

    
    #  dM * phis.
    return  num / den
    
def lumfn_stepwise(vmax, Mcol='MALL_0P0', Mmin_col='DDPMALL_0P0_VISZ', survey='gama'):
    dM        = 0.1
    phi_Ms    = np.arange(-26., -16., dM) 

    # Flat initial
    # phi_init  = 1.e-2 * np.ones_like(phi_Ms)

    # Better convergence; 
    # HACK: shouldn't assume 'truth'
    phi_init  = named_schechter(phi_Ms, named_type='TMR')

    diff      = 1.e99
    phis      = phi_init 

    iteration = 0
    
    while  (diff > 1.e-6):
        print('Solving for iteration {:d} with diff. {:.6e}'.format(iteration, diff))
        
        new_phis = []
    
        for i, (phi_M, phi) in enumerate(zip(phi_Ms, phis)):
            phi_hat = lumfn_stepwise_eval(vmax, phi_M, phi, phis, phi_Ms, dM, Mcol=Mcol, Mmin_col=Mmin_col, survey=survey)
            new_phis.append(phi_hat)

        new_phis    = np.array(new_phis)
        
        diff        = np.sum((new_phis - phis)**2.)
    
        #  Update previous estimate. 
        phis        = new_phis

        iteration  += 1

    phis  = phis / dM 
    nn    = dM * np.sum(phis)
    phis /= nn
    
    return  phi_Ms, phis
        

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Generate Gold stepwise luminosity function.')
    parser.add_argument('-s', '--survey', help='Select survey', default='gama')
    parser.add_argument('--nooverwrite',  help='Do not overwrite outputs if on disk', action='store_true')

    args   = parser.parse_args()

    survey = args.survey
    
    fpath  = findfile('ddp')
    opath  = findfile('lumfn_step')

    print(fpath)
    print(opath)

    if args.nooverwrite:
        overwrite_check(opath)

    ddp    = Table.read(fpath)
    ddp.pprint()

    phi_Ms, phis = lumfn_stepwise(ddp, survey=survey)
    result       = Table(np.c_[phi_Ms, phis], names=['Ms', 'PHI_STEP'])

    result.write(opath, format='fits', overwrite=True)

    print('\n\nDone.\n\n')
