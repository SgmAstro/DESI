from astropy.table import Table

def renormalise_d8LF(cat, fdelta, fdelta_zeropoint, self_count):
    '''
    fscale equal to Equation 7 in McNaught-Roberts (2014).
    See: https://arxiv.org/pdf/1409.4681.pdf
    '''
    
    cat = Table(cat, copy=True)

    if self_count = False:
        cat['PHI_N']           /= fdelta
        cat['PHI_N_ERROR']     /= fdelta
        cat['PHI_IVMAX']       /= fdelta
        cat['PHI_IVMAX_ERROR'] /= fdelta

    else:
        # ADD SELF_COUNT WORK HERE
        
    
         cat['PHI_IVMAX_ERROR'] *= fdelta/fdelta_zeropoint
    ## TODO
    ## DDP1 < M < DDP1
    ## cat['PHI_IVMAX_ERROR'] *= fdelta/fdelta_ddp 

    return cat
