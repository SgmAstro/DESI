import numpy as np

from   astropy.table import Table
from   cosmo import volcom

def vmaxer(dat, zmin, zmax, area, zcol='ZGAMA', extra_cols=[]):
    assert dat[zcol].min() <= zmin
    assert dat[zcol].max() >= zmax

    cols   = [zcol, 'ZMIN', 'ZMAX'] + extra_cols

    result = Table(dat[cols], copy=True)
    result = result[result[zcol] >= zmin]
    result = result[result[zcol] <= zmax]

    result['ZMIN']  = np.clip(result['ZMIN'], zmin, None)
    result['ZMAX']  = np.clip(result['ZMAX'], None, zmax)
    
    result['VMAX']  = volcom(result['ZMAX'], area)
    result['VMAX'] -= volcom(result['ZMIN'], area)

    result['VZ']    = volcom(result[zcol], area)
    result['VZ']   -= volcom(result['ZMIN'], area)
    
    return result
