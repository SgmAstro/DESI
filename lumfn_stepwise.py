import  os
import  sys
import  time
import  tqdm
import  argparse
import  numpy           as     np

from    runtime         import calc_runtime
from    functools       import partial
from    multiprocessing import Pool
from    astropy.table   import Table
from    findfile        import findfile, overwrite_check, write_desitable
from    schechter       import named_schechter
from    ddp             import initialise_ddplimits
from    ddp_limits      import limiting_curve_path


def lum_binner(x, dM):
    '''
    Eqn. 2.10a, W(x), of Efstathiou, Ellis & Peterson.
    '''
    return  np.abs(x) <= (dM / 2.)

def lum_visible(x, dM):
    '''                                                                                                                          
    Eqn. 2.10b, H(x), of Efstathiou, Ellis & Peterson.                                                                            
    
    Note:  not currently used.
    '''
    result    = -x/dM + 1./2.

    result[x >=  (dM / 2.)] = 0.0
    result[x <= -(dM / 2.)] = 1.0
    
    return  result

def process_one(split, Mmins, Mmaxs, dM, phi_Ms, phis):
    '''
    Stepwise weights for each split.
    '''

    weights     = []

    Mmins       = np.array(Mmins[split], copy=True)
    Mmaxs       = np.array(Mmaxs[split], copy=True)

    for i in np.arange(len(Mmins)):
        Mmin    = Mmins[i]
        Mmax    = Mmaxs[i]
            
        isin    = (phi_Ms >= Mmin) & (phi_Ms <= Mmax)

        assert  np.count_nonzero(isin)

        # 1 / <n> weight.
        weight  = 1. / np.sum(phis[isin])
        weight /= dM

        weights.append(weight)

        # print(Mmin, Mmax, weight)

    weights = np.array(weights)

    return  weights.tolist()

def lumfn_stepwise_eval(vmax, dM, phi_M, phi, phi_Ms, phis, Mcol='MCOLOR_0P0', survey='gama', nproc=12):
    '''
    Eqn. 2.12, of Efstathiou, Ellis & Peterson.   
    '''
    vmax        = Table(vmax, copy=True)

    # Fortran indexing, 1 .. 7 inclusive.
    vmax.sort('REST_GMR_0P1_INDEX')

    # MCOLOR_0P0
    Ms          = vmax[Mcol]
    uidxs       = np.unique(vmax['REST_GMR_0P1_INDEX'].data)

    # For each galaxy, namely rest frame gmr_0p1 and z, limiting MCOLOR_0P0 at the bright and faint ends. 
    Mmins       = vmax['STEPWISE_BRIGHTLIM_0P0']
    Mmaxs       = vmax['STEPWISE_FAINTLIM_0P0']
    
    splits      = []
    split_idx   = np.arange(len(vmax))

    nums        = []
    
    # Unbalanced, poor performance for pool. 
    for uidx in uidxs:
        # Solve for this color set. 
        isin    = (vmax['REST_GMR_0P1_INDEX'].data == uidx)        
        isin   &= (lum_binner(vmax[Mcol] - phi_M, dM))

        num     = np.count_nonzero(isin)
        nums.append(num)

        if num == 0:
            splits.append([])
        
        else:            
            # print('{:d}\t{:d}\t{:.6f}\t{:.6f}'.format(uidx, num, vmax[Mcol][isin].min(), phi_M, vmax[Mcol][isin].max()))

            # Volume limited sample for mag. phi_M and this rest-frame color. 
            sub     =  vmax[isin]

            # Liberal limits. 
            zmin    =  sub['ZSURV'].min()
            zmax    =  sub['ZSURV'].max()

            # print(zmin, zmax)

            isin    = (vmax['ZSURV'].data >= zmin) & (vmax['ZSURV'].data <= zmax)

            # print('\t{}\t{}\t{}'.format(zmin, zmax, np.count_nonzero(isin)))
        
            splits.append(split_idx[isin])
    
    results = []

    with Pool(nproc) as pool:
        # for result in tqdm.tqdm(pool.imap(partial(process_one, Mmins=Mmins, Mmaxs=Mmaxs, dM, phi_Ms=phi_Ms, phis=phis), iterable=splits), total=len(splits)):
        for result in pool.imap(partial(process_one, Mmins=Mmins, Mmaxs=Mmaxs, dM=dM, phi_Ms=phi_Ms, phis=phis), iterable=splits):
            results.append(result)
    
        pool.close()

        # https://stackoverflow.com/questions/38271547/when-should-we-call-multiprocessing-pool-join                                                                                                       
        pool.join()

    '''
    for split, result in zip(splits, results):
        if len(split) > 0:
            sub          = vmax[split] 
            sub['INBAR'] = result

            sub.sort('ZSURV')
            sub['ZSURV', 'REST_GMR_0P1', 'REST_GMR_0P1_INDEX', 'MCOLOR_0P0', 'STEPWISE_BRIGHTLIM_0P0', 'STEPWISE_FAINTLIM_0P0', 'INBAR'].pprint()
    '''

    results  = [np.sum(x) for x in results]
    results  =  np.array(results) # [1/dM]

    #  dM * phis?
    nums     = np.array(nums)
    results *= nums
    
    phi_hat  = np.sum(results)

    # print('{:.6f}\t{:.6f}\t{:.6f}\t{:.6f}'.format(phi_M, phi, np.sum(nums), phi_hat))

    return phi_hat

def lumfn_stepwise(vmax, Mcol='MCOLOR_0P0', tolerance=1.e-3):
    # Note: match lumfn binning.
    nbins      = 12  # 36

    phi_Ms     = np.linspace(-23.,  -16.,  12)[::-1]
    dM         = np.abs(np.diff(phi_Ms)[0])

    # Initialise phi estimates - uniform. 
    # phi_init = dM * 1. * np.ones(len(phi_Ms), dtype=float)
    phi_init   = named_schechter(phi_Ms, named_type='TMR')
    
    diff       = 1.e99
    phis       = phi_init 

    iteration  = 0

    # Remove anything not in the limits, as digitize returns 0, len(array) outwith. 
    isin       = (vmax[Mcol] >= phi_Ms.min()) & (vmax[Mcol] <= phi_Ms.max())
    vmax       =  vmax[isin]

    norm       = np.sum(phis)
    
    while (diff > tolerance):
        print('\n\n------------  Solving for iteration {:d} with diff. {:.6e}  ------------'.format(iteration, diff))
        
        new_phis = []
    
        for i, (phi_M, phi) in enumerate(zip(phi_Ms, phis)):
            phi_hat = lumfn_stepwise_eval(vmax, dM, phi_M, phi, phi_Ms, phis, Mcol=Mcol)
            new_phis.append(phi_hat)

        new_phis    = np.array(new_phis)
    
        #  Update previous estimate. 
        _phis       = norm * (new_phis / np.sum(new_phis))

        for phi_M, phi, _phi in zip(phi_Ms, phis, _phis):
            print('{:.6f}\t{:.6e}\t{:.6e}'.format(phi_M, phi, _phi))

        diff        = np.sum((_phis - phis)**2.)
        phis        = _phis

        iteration  += 1

    # phis  = phis / dM 

    # nn    = dM * np.sum(phis)
    # phis /= nn
    
    # print('Final M={} recovers weights for all galaxies in vmax ({} weights for {} galaxies).'.format(phi_M, len(weights), len(vmax)))
    
    return  phi_Ms + dM/2., phis
        

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Generate Gold stepwise luminosity function.')
    parser.add_argument('--log',          help='Create a log file of stdout.', action='store_true')
    parser.add_argument('-s', '--survey', help='Select survey', default='gama')
    parser.add_argument('--dryrun',       help='Dryrun', action='store_true')
    parser.add_argument('--nooverwrite',  help='Do not overwrite outputs if on disk', action='store_true')
    parser.add_argument('--version',      help='Add version', default='GAMA4')
    
    start       = time.time() 

    args        = parser.parse_args()
    log         = args.log
    survey      = args.survey
    dryrun      = args.dryrun
    nooverwrite = args.nooverwrite
    version     = args.version

    if log:
        logfile = findfile(ftype='lumfn_step', dryrun=False, survey=survey, log=True)

        print(f'Logging to {logfile}')

        sys.stdout = open(logfile, 'w')

    fpath       = findfile('ddp',        dryrun=dryrun, survey=survey, version=version)
    opath       = findfile('lumfn_step', dryrun=dryrun, survey=survey, version=version)

    if nooverwrite:
        overwrite_check(opath)

    ddp         = Table.read(fpath)
    ddp.pprint()

    phi_Ms, phis           = lumfn_stepwise(ddp)

    result                 = Table(np.c_[phi_Ms, phis], names=['MID_M', 'PHI_STEPWISE'])

    runtime                = calc_runtime(start, 'Writing {}'.format(opath))    

    result.meta['IMMUTABLE'] = 'FALSE'

    print(f'Writing {opath}')
    
    write_desitable(opath, result)

    runtime = calc_runtime(start, 'Finished')

    if log:
        sys.stdout.close()

