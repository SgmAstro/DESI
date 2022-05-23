import numpy             as     np
import runtime

from   delta8_limits     import d8_limits
from   params            import fillfactor_threshold
from   scipy.interpolate import interp1d
from   cosmo             import volcom
from   ddp_zlimits       import ddp_zlimits
from   astropy.table     import Table, vstack
from   findfile          import findfile, fetch_fields, gather_cat


def volavg_fillfactor(survey='gama', ftype='randoms_bd_ddp_n8', dryrun=False, prefix='randoms_ddp1', write=False, field='G9', tier=None, pprint=False):
    print('\n\nSolving for volume average fillfactor.')

    fields      = fetch_fields(survey=survey)
    rpaths      = [findfile(ftype=ftype, dryrun=dryrun, field=ff, survey=survey, prefix=prefix) for ff in fields]    
    rand        = gather_cat(rpaths)
    nrand       = len(rand)

    if tier != None:
        rand    = rand[rand['DDP1_DELTA8_TIER'].data == tier]

    ## 
    idx         = np.argsort(rand['Z']) 
    sorted_rand = Table(rand, copy=True)[idx]

    del rand

    print('Randoms {:.6f} <= z <= {:.6f}'.format(sorted_rand['Z'].min(), sorted_rand['Z'].max()))

    dbin        = 1.e-3

    zlo         = ddp_zlimits['DDP1'][0]
    zhi         = ddp_zlimits['DDP1'][1]

    bins        = np.arange(zlo, zhi + dbin, dbin)
    idx         = np.digitize(sorted_rand['Z'], bins=bins)

    result      = []

    for i, bb in enumerate(bins):
        sub     = idx <= i
        vfrac   = 1. * np.count_nonzero(sub) / nrand

        sub    &= (sorted_rand['FILLFACTOR'] > fillfactor_threshold)
        cfrac   = 1. * np.count_nonzero(sub) / nrand

        midb    = bb + dbin/2.

        result.append([midb, vfrac, cfrac])

    del sorted_rand

    result      = np.array(result)
    result      = Table(result, names=['Z', 'RANDFRAC', 'RANDFRAC_FILLFACTOR'])
    
    if pprint:
        result.pprint()

    if write:
        opath   = findfile(ftype='volavg_fillfactor', dryrun=dryrun, field=field, survey=survey, prefix=prefix, utier=tier)
        result.write(opath, format='fits', overwrite=True)

    vol_splint  = interp1d(result['Z'], result['RANDFRAC'], kind='linear', copy=True, bounds_error=False, fill_value=0.0, assume_sorted=False)
    cut_splint  = interp1d(result['Z'], result['RANDFRAC_FILLFACTOR'], kind='linear', copy=True, bounds_error=False, fill_value=0.0, assume_sorted=False)

    del result

    return  vol_splint, cut_splint

def eval_volavg_fillfactor(dat, survey='gama', ftype='randoms_bd_ddp_n8', dryrun=False, prefix='randoms_ddp1', write=False, field='G9', tier=None):
    vol_splint, cut_splint = volavg_fillfactor(survey=survey, ftype=ftype, dryrun=dryrun, prefix=prefix, write=write, field=field, tier=tier)

    def _eval_volavg_fillfactor(zmax, zmin):
        return (cut_splint(zmax) - cut_splint(zmin)) / (vol_splint(zmax) - vol_splint(zmin))

    result = []

    for row in dat:
        result.append(_eval_volavg_fillfactor(row['ZMAX'], row['ZMIN']))
        
    result = np.array(result)

    dat['FILLFACTOR_VMAX'] = result

    return  0

def volfracs(rand, bitmasks=['IN_D8LUMFN']):
    '''
    Calculate volume fractions, typically used to rescale VMAX from solid angle and 
    (DDP1) redshift defined to account for e.g. density restrictions and sphere completeness
    restrictions.

    Note:
        At the minute, IN_d8LUMFN is restricted to a fillfactor cut.  As the correction to VMAX
        due to this cut is zmax dependent, we correct the individual galaxy zmaxs for this.  Therefore,                                                                                             
        The volume correction for the randoms should be restricted to a strict d8 volume, having applied
        the fillfactor cut to the randoms already.  Otherwise, we double count this effect (approximately).
    '''
    utiers    = np.unique(rand['DDP1_DELTA8_TIER'].data)
    utiers_zp = np.unique(rand['DDP1_DELTA8_TIER_ZEROPOINT'].data)

    utiers    = utiers[utiers >= 0]
    utiers_zp = utiers_zp[utiers_zp >= 0]

    print('Unique tiers: {}'.format(utiers))

    # Limit randoms to DDP1 redshift limits. 
    ddp1_rand = rand[rand['DDPZLIMS'][:,0] == 1]

    print('DDP1 randoms: {:.6f} < z < {:.6f}'.format(ddp1_rand['Z'].min(), ddp1_rand['Z'].max()))

    for idx in range(3):
        ddp_idx = idx + 1

        # Within a given DDP z limits.
        sub     = rand[rand['DDPZLIMS'][:,idx] == 1]

        rand.meta['DDP{}_FULL8FRAC'.format(ddp_idx)] = np.mean(sub['FILLFACTOR'] > fillfactor_threshold) 
    
    for bm in bitmasks:
        isin      = (ddp1_rand[bm].data == 0)

        # At least a cut on sphere completeness.                                                                                                                                                        
        ddp1_rand = ddp1_rand[isin]

        print(bm, np.mean(isin))

    for ut in range(len(d8_limits)):
        print()

        in_tier = (ddp1_rand['DDP1_DELTA8_TIER'].data == ut)

        # print(ut, d8_limits[ut], np.mean(d8_limits[ut]))

        if np.count_nonzero(in_tier) > 0: 
            rand.meta['DDP1_d{}_VOLFRAC'.format(ut)]   = '{:.6f}'.format(np.mean(in_tier))
            rand.meta['DDP1_d{}_TIERMEDd8'.format(ut)] = '{:.6f}'.format(np.mean(ddp1_rand['DDP1_DELTA8'].data[in_tier]))

        else:
            rand.meta['DDP1_d{}_VOLFRAC'.format(ut)]   = '{:.6f}'.format(0.0)
            rand.meta['DDP1_d{}_TIERMEDd8'.format(ut)] = '{:.6f}'.format(np.mean(d8_limits[ut]))
            
        print('DDP1_d{}_VOLFRAC OF {} added.'.format(ut,    rand.meta['DDP1_d{}_VOLFRAC'.format(ut)]))
        print('DDP1_d{}_TIERMED d8 OF {} added.'.format(ut, rand.meta['DDP1_d{}_TIERMEDd8'.format(ut)]))

        # Zero point.                                                                                                                                                                                      
        in_tier = (ddp1_rand['DDP1_DELTA8_TIER_ZEROPOINT'].data == ut) 

        if np.count_nonzero(in_tier) > 0:
            rand.meta['DDP1_d{}_ZEROPOINT_VOLFRAC'.format(ut)]   = '{:.6f}'.format(np.mean(in_tier))
            rand.meta['DDP1_d{}_ZEROPOINT_TIERMEDd8'.format(ut)] = '{:.6}'.format(np.mean(ddp1_rand['DDP1_DELTA8_ZEROPOINT'].data[in_tier]))
        
        else:
            rand.meta['DDP1_d{}_ZEROPOINT_VOLFRAC'.format(ut)]   = '{:.6f}'.format(0.0)
            rand.meta['DDP1_d{}_ZEROPOINT_TIERMEDd8'.format(ut)] = '{:.6f}'.format(np.mean(d8_limits[ut]))

        print('DDP1_d{}_ZEROPOINT_VOLFRAC OF {:.10f} added.'.format(ut, np.mean(in_tier)))
        print('DDP1_d{}_ZEROPOINT_TIERMED d8 OF {} added.'.format(ut, rand.meta['DDP1_d{}_ZEROPOINT_TIERMEDd8'.format(ut)]))

    return  rand


if __name__ == '__main__':
    fpath = findfile(ftype='randoms_bd_ddp_n8', dryrun=False, field='GALL', survey='gama', prefix='randoms_ddp1')

    rand  = Table.read('/cosma5/data/durham/dc-wils7/GAMA4/randoms/randoms_ddp1_bd_ddp_n8_GALL_0.fits')
    rand.pprint()

    volavg_fillfactor(rand)
