import numpy           as     np

from   astropy.table   import Table, vstack
from   cosmo           import volcom
from   bitmask         import lumfn_mask, consv_mask, update_bit
from   volfracs        import volfracs
from   findfile        import findfile, fetch_fields, write_desitable
from   params          import fillfactor_threshold
from   volfracs        import eval_volavg_fillfactor


def vmaxer_rand(survey='gama', ftype='randoms_bd_ddp_n8', dryrun=False, prefix='randoms_ddp1', bitmasks=['IN_D8LUMFN'], conservative=False, write=False):    
    fields = fetch_fields(survey=survey)

    rpaths = [findfile(ftype=ftype, dryrun=dryrun, field=ff, survey=survey, prefix=prefix) for ff in fields]
    rand   = vstack([Table.read(xx) for xx in rpaths])

    update_bit(rand['IN_D8LUMFN'], lumfn_mask, 'FILLFACTOR', rand['FILLFACTOR'].data < fillfactor_threshold)
    
    rand   = volfracs(rand, bitmasks=bitmasks)    

    # HACK SURVEYHACK
    assert survey == 'gama'

    print('\n\n')

    for key in rand.meta.keys():
        print(key, rand.meta[key])

    print('\n\n')

    if write:
        opath  = findfile(ftype=ftype, dryrun=dryrun, field='GALL', survey=survey, prefix=prefix)

        print(f'Writing {opath}')

        write_desitable(opath, rand)

    return  rand

def vmaxer(dat, zmin, zmax, extra_cols=[], fillfactor=True, bitmasks=['IN_D8LUMFN'], tier=None):
    assert  dat['ZSURV'].min() <= zmin
    assert  dat['ZSURV'].max() >= zmax

    # Columns to be propagated
    extra_cols += ['FIELD', 'MALL_0P0', 'MCOLOR_0P0', 'FIELD', 'IN_D8LUMFN', 'RA', 'DEC', 'DDPMALL_0P0']
    extra_cols += ['FILLFACTOR', 'REST_GMR_0P1_INDEX']


    if 'STEPWISE_BRIGHTLIM_0P0' in dat.dtype.names:
        extra_cols += ['STEPWISE_BRIGHTLIM_0P0', 'STEPWISE_FAINTLIM_0P0']


    if 'WEIGHT_STEPWISE' in dat.dtype.names:
        extra_cols += ['WEIGHT_STEPWISE']
                
    cols        = ['ZSURV', 'ZMIN', 'ZMAX'] + extra_cols
    cols        = list(set(cols))

    print('Solving for VMAX catalog.')

    result      = Table(dat[cols], copy=True)
    result.meta = dat.meta

    # Apply redshift limits.
    result      = result[result['ZSURV'] >= zmin]
    result      = result[result['ZSURV'] <= zmax]

    print('VMAX catalog cut to {:.6f} <= z <= {:.6f}.'.format(zmin, zmax))

    # First, define vmaxes on vl1 redshift limits. 
    result['ZMIN']  = np.clip(result['ZMIN'], zmin, None)
    result['ZMAX']  = np.clip(result['ZMAX'], None, zmax)
    
    area            = dat.meta['AREA']

    result['VZ']    = volcom(result['ZSURV'], area)
    result['VZ']   -= volcom(result['ZMIN'],  area)

    result['VMAX']  = volcom(result['ZMAX'],  area)
    result['VMAX'] -= volcom(result['ZMIN'],  area)

    # Apply bitmask cut. 
    for bmask in bitmasks:
        isin    = result[bmask] == 0
        result  = result[isin]

        print(bmask, np.mean(isin))

    # New limits of subset. 
    zmin        = result['ZSURV'].min()
    zmax        = result['ZSURV'].max()
    
    VV          = volcom(zmax, area) - volcom(zmin, area)

    print('Retrieved area {:.4f} [sq. deg.]'.format(area))

    # Include impact of e.g. bitmask cuts on these quantities.
    result.meta.update({'FORCE_ZMIN': zmin,\
                        'FORCE_ZMAX': zmax,\
                        'FORCE_VOL':    VV})

    # Was fillfactor applied?
    result.meta['FILLFACTOR']     = fillfactor

    if fillfactor:
        # HACK SURVEYHACK
        result['FILLFACTOR_VMAX']  = eval_volavg_fillfactor(result,\
                                                            survey='gama',\
                                                            ftype='randoms_bd_ddp_n8',\
                                                            dryrun=False,\
                                                            prefix='randoms_ddp1',\
                                                            write=False,\
                                                            tier=tier)

        result['VZ']   *= result['FILLFACTOR_VMAX']
        result['VMAX'] *= result['FILLFACTOR_VMAX']
    
    return  result
