import  os
import  time
import  tqdm
import  argparse
import  numpy           as     np

from    runtime         import calc_runtime
from    functools       import partial
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

def process_one(split, Mmins, Mmaxs, phi_Ms, phis):
    weights    = []

    Mmins      = np.array(Mmins[split], copy=True)
    Mmaxs      = np.array(Mmaxs[split], copy=True)

    dM         = np.diff(phi_Ms)[0]

    for i in np.arange(len(Mmins)):
        Mmin   = Mmins[i]
        Mmax   = Mmaxs[i]
            
        isin   = (phi_Ms - (dM / 2.) > Mmin) & (phi_Ms + (dM / 2.) < Mmax)
        weight = 1. / np.sum(phis[isin])

        weights.append(weight)

    weights = np.array(weights) # [dM]

    del Mmins
    del Mmaxs

    return  weights.tolist()

def lumfn_stepwise_eval(vmax, phi_M, phi, phis, phi_Ms, dM, Mcol='MALL_0P0', Mmin_col='DDPMALL_0P0_VISZ', survey='gama', nproc=12):
    '''
    Eqn. 2.12, of Efstathiou, Ellis & Peterson.   
    '''
    Ms    = vmax[Mcol]
    Mmins = vmax[Mmin_col]

    num   = np.count_nonzero(lum_binner(phi_M - Ms))

    bright_curve, bright_curve_r, faint_curve, faint_curve_r = initialise_ddplimits(survey=survey)

    # HACK MALL, MQALL?
    assert  Mcol == 'MALL_0P0' 

    zmin  = bright_curve(phi_M - dM/2.) 

    # HACK: no bright limit; 
    zmin  = 0.0

    zmax  =  faint_curve(phi_M + dM/2.)

    # Deprecated:
    # facs              = lum_visible(phi_M - Mmins) 
    # facs[facs > 0.0] /= np.sum(dM * phis[:,None] * lum_visible(phi_Ms[:,None] - Mmins[facs > 0.0]), axis=0)
    # den               = np.sum(facs) 

    # print(phi_M, zmin, zmax)

    zcol      = 'Z{}'.format(survey.upper())
    vol_lim   = vmax[(vmax[zcol] > zmin) * (vmax[zcol] < zmax)]

    Mmins     = bright_curve_r(vol_lim[zcol])
    
    # HACK: no bright limit. 
    Mmins     = -99. * np.ones_like(Mmins)

    Mmaxs     =  faint_curve_r(vol_lim[zcol])

    if len(vol_lim) == 0:
        return  phi * np.one_like(vol_lim) / len(vol_lim)

    print('{:.4f}\t{:.3f}\t{:.3f}\t{:.2f}\t{:.2f}\t{:d}\t{:d}'.format(phi_M, zmin, zmax, Mmins.min(), Mmaxs.max(), len(vol_lim), len(vmax)))
    
    split_idx = np.arange(len(vol_lim))
    splits    = np.array_split(split_idx, 10 * nproc)

    results   = []

    with Pool(nproc) as pool:
        '''
        result   = pool.imap(partial(process_one, Mmins=Mmins, Mmaxs=Mmaxs, phi_Ms=phi_Ms, phis=phis), iterable=splits)
        results += result
        '''
        
        for result in tqdm.tqdm(pool.imap(partial(process_one, Mmins=Mmins, Mmaxs=Mmaxs, phi_Ms=phi_Ms, phis=phis), iterable=splits), total=len(splits)):
            results += result
    
    results  = np.array(results) # [dM]
    
    # print(result)

    results *= dM

    pool.close()

    # print(len(results), len(vol_lim))

    #  Deprecated
    #  return  num / den

    #  dM * phis.   
    return  results

def lumfn_stepwise(vmax, Mcol='MALL_0P0', Mmin_col='DDPMALL_0P0_VISZ', survey='gama', tolerance=1.):
    dM        = 0.5
    phi_Ms    = np.arange(-25.5, -15.5, dM)[::-1]

    print(f'Solving for {phi_Ms}')

    print('{}\t{}'.format(phi_Ms.max(), vmax[Mcol].max()))
    print('{}\t{}'.format(phi_Ms.min(), vmax[Mcol].min()))
    
    # Flat initial
    # phi_init  = 1.e-2 * np.ones_like(phi_Ms)

    # Better convergence; 
    # HACK: shouldn't assume 'truth'
    phi_init  = named_schechter(phi_Ms, named_type='TMR')

    diff      = 1.e99
    phis      = phi_init 

    iteration = 0
    
    while  (diff > tolerance):
        print('Solving for iteration {:d} with diff. {:.6e}'.format(iteration, diff))
        
        new_phis    = []
    
        for i, (phi_M, phi) in enumerate(zip(phi_Ms, phis)):
            weights = lumfn_stepwise_eval(vmax, phi_M, phi, phis, phi_Ms, dM, Mcol=Mcol, Mmin_col=Mmin_col, survey=survey)
            
            phi_hat = np.sum(weights)
            new_phis.append(phi_hat)

        new_phis    = np.array(new_phis)
        
        diff        = np.sum((new_phis - phis)**2.)

        #  HACK
        diff        = 0.1 * tolerance
    
        #  Update previous estimate. 
        phis        = new_phis

        iteration  += 1

    phis  = phis / dM 
    nn    = dM * np.sum(phis)
    phis /= nn

    print('Final M={} recovers weights for all galaxies in vmax ({} weights for {} galaxies).'.format(phi_M, len(weights), len(vmax)))
    
    return  phi_Ms, phis, weights
        

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Generate Gold stepwise luminosity function.')
    parser.add_argument('-s', '--survey', help='Select survey', default='gama')
    parser.add_argument('--nooverwrite',  help='Do not overwrite outputs if on disk', action='store_true')

    start  = time.time() 

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

    phi_Ms, phis, weights = lumfn_stepwise(ddp, survey=survey)
    result                = Table(np.c_[phi_Ms, phis], names=['Ms', 'PHI_STEP'])

    runtime = calc_runtime(start, 'Writing {}'.format(opath))    
    result.write(opath, format='fits', overwrite=True)

    ddp                    = Table.read(fpath) 
    ddp['WEIGHT_STEPWISE'] = weights

    ddp.write(fpath, format='fits', overwrite=True)

    runtime = calc_runtime(start, 'Finished')
