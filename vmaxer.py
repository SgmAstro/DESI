import numpy           as     np

from   astropy.table   import Table, vstack
from   cosmo           import volcom
from   bitmask         import lumfn_mask, consv_mask, update_bit
from   volfracs        import volfracs
from   findfile        import findfile, fetch_fields


def vmaxer_rand(survey='gama', ftype='randoms_bd_ddp_n8', dryrun=False, prefix='randoms_ddp1', bitmasks=['IN_D8LUMFN'], conservative=False):    
    fields = fetch_fields(survey=survey)

    rpaths = [findfile(ftype=ftype, dryrun=dryrun, field=ff, survey=survey, prefix=prefix) for ff in fields]
    rand   = vstack([Table.read(xx) for xx in rpaths])

    update_bit(rand['IN_D8LUMFN'], lumfn_mask, 'FILLFACTOR', rand['FILLFACTOR'].data < 0.8)
    
    rand   = volfracs(rand, bitmasks=bitmasks)    

    return  rand

def vmaxer(dat, zmin, zmax, extra_cols=[], fillfactor=True, bitmasks=['IN_D8LUMFN']):
    assert  dat['ZSURV'].min() <= zmin
    assert  dat['ZSURV'].max() >= zmax

    # Columns to be propagated
    extra_cols += ['MALL_0P0', 'MCOLOR_0P0', 'FIELD', 'IN_D8LUMFN', 'RA', 'DEC']

    if 'WEIGHT_STEPWISE' in dat.dtype.names:
        extra_cols += ['WEIGHT_STEPWISE']
        
    if fillfactor == True:
        extra_cols += ['FILLFACTOR', 'FILLFACTOR_VMAX']
        
    cols        = ['ZSURV', 'ZMIN', 'ZMAX'] + extra_cols
    cols        = list(set(cols))

    result      = Table(dat[cols], copy=True)
    result.meta = dat.meta

    # Apply redshift limits.
    result      = result[result['ZSURV'] >= zmin]
    result      = result[result['ZSURV'] <= zmax]
    
    # Apply bitmask cut. 
    for bmask in bitmasks:
        isin    = result[bmask] == 0
        result  = result[isin]

        print(bmask, np.mean(isin))

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

    result.meta['FILLFACTOR']     = fillfactor

    if fillfactor:
        # Wrap volavg fillfactor(< z) required for vmax.
        fillfactor_vmax_min       = result['FILLFACTOR_VMAX'].min()
        fillfactor_vmax_max       = result['FILLFACTOR_VMAX'].max()

        result['FILLFACTOR_VMAX'] = np.clip(result['FILLFACTOR_VMAX'], fillfactor_vmax_min, fillfactor_vmax_max)

        result['VMAX']           *= result['FILLFACTOR_VMAX']
    
    return  result
