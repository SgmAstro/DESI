import  os
import  sys
import  time
import  tqdm
import  argparse
import  numpy           as     np
import  multiprocessing

from    runtime         import calc_runtime
from    functools       import partial
from    multiprocessing import Pool
from    astropy.table   import Table
from    findfile        import findfile, overwrite_check, write_desitable
from    schechter       import named_schechter
from    ddp             import initialise_ddplimits
from    ddp_limits      import limiting_curve_path
from    ddp_zlimits     import ddp_zlimits
from    params          import fillfactor_threshold
from    schechter       import named_schechter


def lum_binner(x, dM):
    '''
    Eqn. 2.10a, W(x), of Efstathiou, Ellis & Peterson.
    '''
    return  np.abs(x) <= (dM / 2.)

def lum_visible(x, dM):
    '''                                                                                                                          
    Eqn. 2.10b, H(x), of Efstathiou, Ellis & Peterson.                                                                            
    
    Note:  
        Not currently used.
    '''
    result    = -x/dM + 1./2.

    result[x >=  (dM / 2.)] = 0.0
    result[x <= -(dM / 2.)] = 1.0
    
    return  result

def process_one(split, Mmins, Mmaxs, dM, phi_Ms, phis):
    '''
    Stepwise (1 / <n> @ z, rest gmr, etc.) weights for each split.
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
        nbar    = dM * np.sum(phis[isin])
        weight  = 1. / nbar

        weights.append(weight)

        # print(Mmin, Mmax, weight)

    weights = np.array(weights)

    return  weights.tolist()

def lumfn_stepwise_eval(vmax, dM, phi_M, phi, phi_Ms, phis, Mcol='MCOLOR_0P0', survey='gama', nproc=12):
    '''
    Eqn. 2.12, of Efstathiou, Ellis & Peterson.   
    '''

    print(f'\n\n----------  Solving for phi_M {phi_M}  ----------')

    # Fortran indexing, 1 .. 7 inclusive.
    # vmax.sort('REST_GMR_0P1_INDEX')

    # HACK:  Color independent
    # vmax['REST_GMR_0P1_INDEX'] = 1

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
            Mmin    =  vmax[Mcol][isin].data.min()
            Mmax    =  vmax[Mcol][isin].data.max()

            # zmin  =  vmax['ZSURV'][isin].data.min()
            # zmax  =  vmax['ZSURV'][isin].data.max()

            zmin    =  vmax['ZMIN'][isin].data.min()
            zmax    =  vmax['ZMAX'][isin].data.max()
            
            print('{}\t{}\t{}\t{}\t{}\t{}'.format(phi_M, Mmin, Mmax, zmin, zmax, np.count_nonzero(isin)))

            isin    = (vmax['ZSURV'].data >= zmin) & (vmax['ZSURV'].data <= zmax)
            
            splits.append(split_idx[isin])
    
    results = []

    with multiprocessing.get_context('spawn').Pool(nproc) as pool:
        # For this phi_M, per rest frame color list of the stepwise (1/<n>) weight for all galaxies in the vol. limited sample.
        for result in pool.imap(partial(process_one, Mmins=Mmins, Mmaxs=Mmaxs, dM=dM, phi_Ms=phi_Ms, phis=phis), iterable=splits):
            results.append(np.array(result))
    
        pool.close()

        # https://stackoverflow.com/questions/38271547/when-should-we-call-multiprocessing-pool-join                                                                                                       
        pool.join()

    '''
    for split, result in zip(splits, results):
        if len(split) > 0:
            sub         = vmax[split] 
            sub['NBAR'] = 1. / result

            sub.sort('ZSURV')
            sub['ZSURV', 'REST_GMR_0P1', 'REST_GMR_0P1_INDEX', 'MCOLOR_0P0', 'STEPWISE_BRIGHTLIM_0P0', 'STEPWISE_FAINTLIM_0P0', 'NBAR'].pprint()
    '''

    results             = np.array([np.sum(x) for x in results])

    nums                = 1. * np.array(nums)
    results[nums == 0.] = 1.

    # dM * phis 
    phi_hat_nocolor     = np.sum(nums) / np.sum(results) 

    # dM * phis 
    phi_hat             = np.sum(nums / results)

    print('{:.6f}\t{:.6e}\t{:.6f}\t{:.6e}\t{:.6e}\t{:.6e}'.format(phi_M, phi, np.sum(nums), np.sum(results), phi_hat_nocolor, phi_hat))

    return  phi_hat, np.sum(nums)

def lumfn_stepwise(vmax, Mcol='MCOLOR_0P0', tolerance=1.e-3, d8=None, normalise=True):
    # Note: match lumfn binning.
    nbins      = 36

    phi_Ms     = np.linspace(-23.,  -16.,  nbins)
    dM         = np.abs(np.diff(phi_Ms)[0])

    # Initialise phi estimates - uniform. 
    phi_inits  = dM * named_schechter(phi_Ms + dM/2., named_type='TMR')
    
    phis       = phi_inits
    norm       = np.sum(phis)

    if d8 != None:
        norm  *= (1. + d8) / (1. + 0.007)

    iteration  = 0
    diff       = 1.e99

    # Remove anything not in the limits, as digitize returns 0, len(array) outwith. 
    isin       = (vmax[Mcol] >= phi_Ms.min()) & (vmax[Mcol] <= phi_Ms.max())
    vmax       =  vmax[isin]

    while (diff > tolerance):
        print('\n\n------------  Solving for iteration {:d} with diff. {:.6e}  ------------'.format(iteration, diff))

        nMs      = []
        new_phis = []
    
        for i, (phi_M, phi) in enumerate(zip(phi_Ms, phis)):
            phi_hat, nM = lumfn_stepwise_eval(vmax, dM, phi_M, phi, phi_Ms, phis, Mcol=Mcol)
            
            nMs.append(nM)
            new_phis.append(phi_hat)

        nMs         = np.array(nMs)
        new_phis    = np.array(new_phis)
    
        #  Update previous estimate. 
        if normalise:
            _phis   = norm * (new_phis / np.sum(new_phis))
        
        else:
            _phis   = new_phis

        print('\n\n------------  Solved for iteration {:d}  ------------'.format(iteration))

        for nM, phi_M, phi_init, phi, _phi in zip(nMs, phi_Ms, phi_inits, phis, _phis):
            print('{:.3f}\t{:.6f}\t{:.6f}\t{:.6f}\t{:.6f}'.format(nM, phi_M, np.log10(phi_init), np.log10(phi), np.log10(_phi)))

        diff        = np.sum((_phis - phis)**2.)
        phis        = _phis

        iteration  += 1

    if normalise:
        isin = (nMs >= 5) & np.isfinite(phis)
        norm = np.sum(phi_inits[isin])

        if d8 != None:
            norm  *= (1. + d8) / (1. + 0.007)

        # print(norm, np.sum(phis[isin]))

        phis *= (norm / np.sum(phis[isin]))

    phis    = phis / dM 
    phi_Ms += dM/2.

    # print('Final M={} recovers weights for all galaxies in vmax ({} weights for {} galaxies).'.format(phi_M, len(weights), len(vmax)))
    
    result_stepwise                       = Table(np.c_[phi_Ms, phis, nMs], names=['MID_M', 'PHI_STEPWISE', 'N'])
    result_stepwise['VALID']              = result_stepwise['N'] >= 5 
    result_stepwise['REF_SCHECHTER']      = named_schechter(result_stepwise['MID_M'], named_type='TMR')

    if d8 != None:
        # TODO HARDCODE 0.007                                                                                                                                                                              
        result_stepwise['REF_SCHECHTER'] *= (1. + d8) / (1. + 0.007)
    
    result_stepwise['REF_RATIO']          = result_stepwise['PHI_STEPWISE'] / result_stepwise['REF_SCHECHTER']

    result_stepwise.meta['DDP1_D8']       = d8
    result_stepwise.meta['EXTNAME']       = 'LUMFN_STEP'
    
    return  result_stepwise
        

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
    
    zlo         = ddp_zlimits['DDP1'][0]
    zhi         = ddp_zlimits['DDP1'][1] 

    isin        = (ddp['ZSURV'] >= zlo) & (ddp['ZSURV'] <= zhi)

    ddp         = ddp[isin]
    ddp.pprint()

    ddp['ZMIN']            = np.clip(ddp['ZMIN'], zlo, None)
    ddp['ZMAX']            = np.clip(ddp['ZMAX'], None, zhi)

    result                 = lumfn_stepwise(ddp)
    '''
    runtime                = calc_runtime(start, 'Writing {}'.format(opath))    

    print(f'Writing {opath}')
    
    write_desitable(opath, result)
    '''
    runtime                = calc_runtime(start, 'Finished')

    if log:
        sys.stdout.close()

