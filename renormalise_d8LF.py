from astropy.table import Table

def renormalise_d8LF(cat, fdelta):
    '''
    fscale equal to Equation 7 in McNaught-Roberts (2014).
    See: https://arxiv.org/pdf/1409.4681.pdf
    '''
    
    cat = Table(cat, copy=True)

    cat['PHI_N']           /= fdelta
    cat['PHI_N_ERROR']     /= fdelta
    cat['PHI_IVMAX']       /= fdelta
    cat['PHI_IVMAX_ERROR'] /= fdelta
    
    return cat
