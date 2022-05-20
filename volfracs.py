import numpy         as     np

from   delta8_limits import d8_limits
from   params        import fillfactor_threshold


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
