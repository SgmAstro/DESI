import numpy as np

from   astropy.table import Table
from   cosmo import volcom


def vmaxer(dat, zmin, zmax, zcol='ZGAMA', extra_cols=[], rand=None):
    assert  dat[zcol].min() <= zmin
    assert  dat[zcol].max() >= zmax

    if rand != None:
        extra_cols += ['FILLFACTOR']

    cols        = [zcol, 'ZMIN', 'ZMAX'] + extra_cols

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

    if rand != None:
        vmax_rand                        = rand[(zmin < rand['Z']) & (rand['Z'] < zmax)]

        result['IN_SAMPLE']              = result['FILLFACTOR'].data > 0.8
        result.meta['IN_SAMPLE_VOLFRAC'] = np.mean(vmax_rand['FILLFACTOR'].data > 0.8)

        print('Fraction of galaxies making fillfactor cut: {:.4f}; fraction of volume (randoms): {:.4f}'.format(np.mean(result['IN_SAMPLE']),\
                                                                                                                result.meta['IN_SAMPLE_VOLFRAC']))

    result['ZMIN']  = np.clip(result['ZMIN'], zmin, None)
    result['ZMAX']  = np.clip(result['ZMAX'], None, zmax)
    
    result['VMAX']  = volcom(result['ZMAX'], area)
    result['VMAX'] -= volcom(result['ZMIN'], area)

    result['VZ']    = volcom(result[zcol], area)
    result['VZ']   -= volcom(result['ZMIN'], area)

    return  result
