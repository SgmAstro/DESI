from astropy.table import Table
from ddp import tmr_DDP1

def renormalise_d8LF(cat, fdelta, fdelta_zeropoint, self_count=False):
    '''
    fscale equal to Equation 7 in McNaught-Roberts (2014).
    See: https://arxiv.org/pdf/1409.4681.pdf
    '''

    print('Renormalising LF with fdelta: {:.6e} and fdelta_zeropoint:  {:.6e}'.format(fdelta, fdelta_zeropoint))
    
    cat = Table(cat, copy=True)
    
    cat['PHI_N']           /= fdelta
    cat['PHI_N_ERROR']     /= fdelta
    cat['PHI_IVMAX']       /= fdelta
    cat['PHI_IVMAX_ERROR'] /= fdelta
    
    if self_count:
        # tmr_DDP1: [-21.8, -20.1] 
        is_ddp1 = (cat['MEDIAN_M'] > tmr_DDP1[0]) & (cat['MEDIAN_M'] < tmr_DDP1[1])

        for col in ['PHI_N', 'PHI_N_ERROR', 'PHI_IVMAX', 'PHI_IVMAX_ERROR']:
            cat[col][is_ddp1] *= (fdelta / fdelta_zeropoint)

    return  cat
