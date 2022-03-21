import numpy         as     np

from   astropy.table import Table
from   cosmo         import volcom

def vmaxer(dat, zmin, zmax, zcol, extra_cols=[], fillfactor=True, conservative=False):
    assert  dat['ZSURV'].min() <= zmin
    assert  dat['ZSURV'].max() >= zmax
        
    extra_cols += ['MCOLOR_0P0', 'FIELD', 'WEIGHT_STEPWISE', 'IN_D8LUMFN']

    if rand is not None:
        extra_cols += ['FILLFACTOR', 'FILLFACTOR_VMAX']

    cols        = ['ZSURV', 'ZMIN', 'ZMAX'] + extra_cols
    cols        = list(set(cols))

    result      = Table(dat[cols], copy=True)
    result.meta = dat.meta

    result      = result[result['ZSURV'] >= zmin]
    result      = result[result['ZSURV'] <= zmax]
    
    # I think this is better incorporated in lumfn
    # ddp.fits does not have CONSERVATIVE
    #if conservative:
    #    result  = result[result['CONSERVATIVE'] == 0] 

    zmin        = result['ZSURV'].min()
    zmax        = result['ZSURV'].max()

    area        = dat.meta['AREA']
    VV          = volcom(zmax, area) - volcom(zmin, area)

    print('Retrieved area {:.4f} [sq. deg.]'.format(area))
    
    result.meta.update({'FORCE_ZMIN': zmin,\
                        'FORCE_ZMAX': zmax,\
                        'VOLUME':       VV})

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
    result.meta['FILLFACTOR'] = fillfactor

    if fillfactor:
        result          = result[result['IN_D8LUMFN'] == 0]
        result['VMAX'] *= result['FILLFACTOR_VMAX']
    
    return  result
