from astropy.table import Table

def renormalise_d8LF(cat, fdelta, fdelta_zeropoint, mags=[-21.8, -20.1], self_count=False):
    '''
    fscale equal to Equation 7 in McNaught-Roberts (2014).
    See: https://arxiv.org/pdf/1409.4681.pdf
    '''
    
    cat = Table(cat, copy=True)

    assert mags[0] < mags[1]
    
    mag_mask = (cat['MEDIAN_M'] > mags[0]) & (result['MEDIAN_M'] < mags[1])
    
    if self_count = True:
        cat[~mag_mask]['PHI_N']           /= fdelta
        cat[~mag_mask]['PHI_N_ERROR']     /= fdelta
        cat[~mag_mask]['PHI_IVMAX']       /= fdelta
        cat[~mag_mask]['PHI_IVMAX_ERROR'] /= fdelta

        cat[mag_mask]['PHI_N']           /= fdelta_zeropoint
        cat[mag_mask]['PHI_N_ERROR']     /= fdelta_zeropoint
        cat[mag_mask]['PHI_IVMAX']       /= fdelta_zeropoint
        cat[mag_mask]['PHI_IVMAX_ERROR'] /= fdelta_zeropoint
        
        
    else:
        cat['PHI_N']           /= fdelta
        cat['PHI_N_ERROR']     /= fdelta
        cat['PHI_IVMAX']       /= fdelta
        cat['PHI_IVMAX_ERROR'] /= fdelta
        
    ## TODO
    ## DDP1 < M < DDP1
    ## cat['PHI_IVMAX_ERROR'] *= fdelta/fdelta_ddp 

    return cat
