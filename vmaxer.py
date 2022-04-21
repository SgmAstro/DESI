import numpy           as     np

from   astropy.table   import Table, vstack
from   cosmo           import volcom
from   bitmask         import lumfn_mask, consv_mask
from   volfracs        import volfracs
from   findfile        import findfile, fetch_fields

def vmaxer_rand(survey='gama', ftype='randoms_bd_ddp_n8', dryrun=False, prefix='', conservative=False, version='GAMA4'):    
    fields = fetch_fields(survey=survey)

    rpaths = [findfile(ftype=ftype, dryrun=dryrun, field=ff, survey=survey, prefix=prefix, version=version) for ff in fields]
    #rand   = [Table.read(xx) for xx in rpaths]

    rand = Table.read(rpaths[0])
    for xx in range(1, len(rpaths)):
        rand_stack = Table.read(rpaths[xx])
        rand = vstack([rand, rand_stack])

    #rand   = rand[rand['ZSURV'] >= zmin]
    #rand   = rand[rand['ZSURV'] <= zmax]

    if conservative:
        rand['CONSERVATIVE'] += (rand['BOUND_DIST'].data < 8.) * consv_mask.BOUNDDIST

        isin                  = (rand['ZSURV'] < 0.9 * rand.meta['DDP1_ZMAX']) & (rand['ZSURV'] > 1.1 * rand.meta['DDP1_ZMIN'])
        rand['CONSERVATIVE'][~isin] += consv_mask.DDP1ZLIM

    # TODO: volfracs currently assumes fillfactor > 0.8 rather than more general bit cut.                                                                                                                 
    rand = volfracs(rand)

    # TODO: define fdelta and d8 based on this all-field rand. 

    
    '''
    # Deprecated: 
    # 
    # 
    # Calculated for DDP1 redshift limits.                                                                                                                                                                
    fdelta = np.array([float(x.meta['DDP1_d{}_VOLFRAC'.format(idx)]) for x in all_rands])
    d8     = np.array([float(x.meta['DDP1_d{}_TIERMEDd8'.format(idx)]) for x in all_rands])

    fdelta_zeropoint = np.array([float(x.meta['DDP1_d{}_ZEROPOINT_VOLFRAC'.format(idx)]) for x in all_rands])
    d8_zeropoint     = np.array([float(x.meta['DDP1_d{}_ZEROPOINT_TIERMEDd8'.format(idx)]) for x in all_rands])

    print('Field vol renormalization: {}'.format(fdelta))
    print('Field d8  renormalization: {}'.format(d8))

    fdelta = fdelta.mean()
    d8     = d8.mean()

    fdelta_zeropoint = fdelta_zeropoint.mean()
    d8_zeropoint     = d8_zeropoint.mean()

    print('Found mean vol. renormalisation scale of {:.3f}'.format(fdelta))
    print('Found mean  d8  renormalisation scale of {:.3f}'.format(d8))
    '''

    # raise  NotImplementedError()

    return  rand

def vmaxer(dat, zmin, zmax, extra_cols=[], fillfactor=True, conservative=False):
    assert  dat['ZSURV'].min() <= zmin
    assert  dat['ZSURV'].max() >= zmax

    # Columns to be propagated
    extra_cols += ['MALL_0P0', 'MCOLOR_0P0', 'FIELD', 'WEIGHT_STEPWISE', 'IN_D8LUMFN']

    if fillfactor == True:
        extra_cols += ['FILLFACTOR', 'FILLFACTOR_VMAX']
        
    if conservative == True:
        extra_cols += ['CONSERVATIVE']


    cols        = ['ZSURV', 'ZMIN', 'ZMAX'] + extra_cols
    cols        = list(set(cols))

    result      = Table(dat[cols], copy=True)
    result.meta = dat.meta

    # Apply redshift limits.
    result      = result[result['ZSURV'] >= zmin]
    result      = result[result['ZSURV'] <= zmax]
    
    result.meta.update({'FORCE_ZMIN': zmin,\
                        'FORCE_ZMAX': zmax})

    # New limits of subset. 
    zmin        = result['ZSURV'].min()
    zmax        = result['ZSURV'].max()

    area        = dat.meta['AREA']
    
    VV          = volcom(zmax, area) - volcom(zmin, area)

    print('Retrieved area {:.4f} [sq. deg.]'.format(area))
    
    result.meta.update({'VOLUME': VV})

    result['ZMIN']  = np.clip(result['ZMIN'], zmin, None)
    result['ZMAX']  = np.clip(result['ZMAX'], None, zmax)
    
    result['VMAX']  = volcom(result['ZMAX'], area)
    result['VMAX'] -= volcom(result['ZMIN'], area)

    result['VZ']    = volcom(result['ZSURV'], area)
    result['VZ']   -= volcom(result['ZMIN'], area)

    result.meta['CONSERVATIVE'] = conservative
    result.meta['FILLFACTOR']   = fillfactor

    if fillfactor:
        # Wrap volavg fillfactor(< z) required for vmax.   
        # TODO:  assumes monotonic.
        fillfactor_vmax_min       = result['FILLFACTOR_VMAX'][result['ZSURV'] >= zmin].min()
        fillfactor_vmax_max       = result['FILLFACTOR_VMAX'][result['ZSURV'] <= zmax].max()

        result['FILLFACTOR_VMAX'] = np.clip(result['FILLFACTOR_VMAX'], fillfactor_vmax_min, fillfactor_vmax_max)

        result['VMAX'] *= result['FILLFACTOR_VMAX']
    
    if conservative:
        result['CONSERVATIVE'] += (result['BOUND_DIST'].data < 8.) * consv_mask.BOUNDDIST

        isin                    = (result['ZSURV'] < 0.9 * result.meta['DDP1_ZMAX']) & (dat['ZSURV'] > 1.1 * result.meta['DDP1_ZMIN'])                                                                    
        result['CONSERVATIVE'][~isin] += consv_mask.DDP1ZLIM 

    return  result
