import numpy           as     np

from   astropy.table   import Table
from   cosmo           import volcom
from   bitmask         import BitMask, lumfn_mask, consv_mask
from   gen_rand_ddp_n8 import volfracs()


def vmaxer(dat, zmin, zmax, zcol, extra_cols=[], fillfactor=True, conservative=False):
    assert  dat['ZSURV'].min() <= zmin
    assert  dat['ZSURV'].max() >= zmax

    # Columns to be propagated
    extra_cols += ['MCOLOR_0P0', 'FIELD', 'WEIGHT_STEPWISE', 'IN_D8LUMFN', 'CONSERVATIVE']

    if rand is not None:
        extra_cols += ['FILLFACTOR', 'FILLFACTOR_VMAX']

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

    # Wrap volavg fillfactor(< z) required for vmax.   
    # TODO:  assumes monotonic.
    fillfactor_vmax_min       = result['FILLFACTOR_VMAX'][result['Z'] >= zmin].min()
    fillfactor_vmax_max       = result['FILLFACTOR_VMAX'][result['Z'] <= zmax].max()
        
    result['FILLFACTOR_VMAX'] = np.clip(result['FILLFACTOR_VMAX'], fillfactor_vmax_min, fillfactor_vmax_max)

    result['ZMIN']  = np.clip(result['ZMIN'], zmin, None)
    result['ZMAX']  = np.clip(result['ZMAX'], None, zmax)
    
    result['VMAX']  = volcom(result['ZMAX'], area)
    result['VMAX'] -= volcom(result['ZMIN'], area)

    result['VZ']    = volcom(result['ZSURV'], area)
    result['VZ']   -= volcom(result['ZMIN'], area)

    result.meta['CONSERVATIVE'] = conservative
    result.meta['FILLFACTOR']   = fillfactor

    if fillfactor:
        result['VMAX'] *= result['FILLFACTOR_VMAX']
    
    if conservative:
        result['CONSERVATIVE'] += (result['BOUND_DIST'].data < 8.) * consv_mask.BOUNDDIST

        isin                    = (result['ZSURV'] < 0.9 * result.meta['DDP1_ZMAX']) & (dat['ZSURV'] > 1.1 * result.meta['DDP1_ZMIN'])                                                                    
        result['CONSERVATIVE'][~isin] += consv_mask.DDP1ZLIM 

    # TODO: volfracs currently assumes fillfactor > 0.8 rather than more general bit cut. 
    result = volfracs(result)
        
    return  result
