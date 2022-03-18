import numpy         as     np

from   astropy.table import Table
from   cosmo         import volcom

def vmaxer(dat, zmin, zmax, zcol='ZGAMA', extra_cols=[], rand=None):
    assert  dat[zcol].min() <= zmin
    assert  dat[zcol].max() >= zmax

    extra_cols += ['MCOLOR_0P0', 'FIELD', 'WEIGHT_STEPWISE', 'IN_D8LUMFN']

    if rand is not None:
        extra_cols += ['FILLFACTOR', 'FILLFACTOR_VMAX']

    cols        = [zcol, 'ZMIN', 'ZMAX'] + extra_cols
    cols        = list(set(cols))

    area        = dat.meta['AREA']
    VV          = volcom(zmax, area) - volcom(zmin, area)

    print('Retrieved area {:.4f} [sq. deg.]'.format(area))
    
    result      = Table(dat[cols], copy=True)

    result.meta = dat.meta
    result.meta.update({'FORCE_ZMIN': zmin,\
                        'FORCE_ZMAX': zmax,\
                        'VOLUME':       VV})

    result      = result[result[zcol] >= zmin]
    result      = result[result[zcol] <= zmax]

    if rand is not None:
        vmax_rand                        = rand[(zmin < rand['Z']) & (rand['Z'] < zmax)]

        fillfactor_vmax_min              = result['FILLFACTOR_VMAX'][result['ZMAX'] >= zmin].min()
        fillfactor_vmax_max              = result['FILLFACTOR_VMAX'][result['ZMAX'] <= zmax].max()
        
        result['FILLFACTOR_VMAX']        = np.clip(result['FILLFACTOR_VMAX'], fillfactor_vmax_min, fillfactor_vmax_max)

    result['ZMIN']  = np.clip(result['ZMIN'], zmin, None)
    result['ZMAX']  = np.clip(result['ZMAX'], None, zmax)
    
    result['VMAX']  = volcom(result['ZMAX'], area)
    result['VMAX'] -= volcom(result['ZMIN'], area)

    result['VZ']    = volcom(result[zcol], area)
    result['VZ']   -= volcom(result['ZMIN'], area)

    return  result
