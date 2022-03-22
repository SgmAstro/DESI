import numpy as np

def volfracs(rand, bitmasks=[]):
    utiers    = np.unique(rand['DDP1_DELTA8_TIER'].data)
    utiers_zp = np.unique(rand['DDP1_DELTA8_TIER_ZEROPOINT'].data)

    print('Unique tiers: {}'.format(utiers))

    ddp1_rand = rand[rand['DDPZLIMS'][:,0]]
    
    for ut in utiers:
        in_tier = (ddp1_rand['DDP1_DELTA8_TIER'].data == ut)

        for bm in bitmasks:
            # Deprecated: (ddp1_rand['FILLFACTOR'].data >= 0.8)                                                                                                                                             
            in_tier &= (ddp1_rand[bm].data == 0)

        rand.meta['DDP1_d{}_VOLFRAC'.format(ut)]   = '{:.6e}'.format(np.mean(in_tier))
        rand.meta['DDP1_d{}_TIERMEDd8'.format(ut)] = '{:.6e}'.format(np.median(ddp1_rand['DDP1_DELTA8'].data[in_tier]))

        print('DDP1_d{}_VOLFRAC OF {:.4f} added.'.format(ut, np.mean(in_tier)))
        print('DDP1_d{}_TIERMED d8 OF {} added.'.format(ut, rand.meta['DDP1_d{}_TIERMEDd8'.format(ut)]))

        # Zero point.                                                                                                                                                                                       
        in_tier = (ddp1_rand['DDP1_DELTA8_TIER_ZEROPOINT'].data == ut) & (ddp1_rand['FILLFACTOR'].data >= 0.8)

        rand.meta['DDP1_d{}_ZEROPOINT_VOLFRAC'.format(ut)]   = '{:.6e}'.format(np.mean(in_tier))
        rand.meta['DDP1_d{}_ZEROPOINT_TIERMEDd8'.format(ut)] = '{:.6e}'.format(np.median(ddp1_rand['DDP1_DELTA8_ZEROPOINT'].data[in_tier]))

        print('DDP1_d{}_ZEROPOINT_VOLFRAC OF {:.4f} added.'.format(ut, np.mean(in_tier)))
        print('DDP1_d{}_ZEROPOINT_TIERMED d8 OF {} added.'.format(ut, rand.meta['DDP1_d{}_ZEROPOINT_TIERMEDd8'.format(ut)]))

    return rand
