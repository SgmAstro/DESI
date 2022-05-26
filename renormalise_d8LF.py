import numpy         as     np

from   astropy.table import Table
from   ddp           import tmr_DDP1


def renormalise_d8LF(idx, lf, fdelta, fdelta_zeropoint, self_count=False, Mcol='MID_M'):
    '''
    fscale equal to Equation 7 in McNaught-Roberts (2014).
    See: https://arxiv.org/pdf/1409.4681.pdf
    '''

    print('Renormalising (IVMAX) LF with log10(fdelta)={:.6e} for tier {}'.format(np.log10(fdelta), idx))
    
    lf    = Table(lf, copy=True)
    
    cols  = ['PHI_N', 'PHI_N_ERROR', 'PHI_IVMAX', 'PHI_IVMAX_ERROR']
    cols += ['PHI_IVMAX_JK', 'PHI_IVMAX_ERROR_JK']
    cols += ['PHI_STEPWISE', 'REF_RATIO']

    cols  = [x for x in cols if x in lf.dtype.names]

    print('Renormalising ...')

    for col in cols:
        print(f'\t{col}')

        lf[col] /= fdelta
    
    if self_count:
        print('Applying log10|self-count correction| of {:.6f}'.format(np.log10(fdelta / fdelta_zeropoint)))
            
        # tmr_DDP1: [-21.8, -20.1] 
        is_ddp1 = (lf[Mcol] > tmr_DDP1[0]) & (lf[Mcol] < tmr_DDP1[1])

        for col in cols:
            lf[col][is_ddp1] *= (fdelta / fdelta_zeropoint)

    lf.meta['SELFCOUNT']       = self_count
    lf.meta['DDP1_VOLFRAC']    = fdelta
    lf.meta['DDP1_VOLFRAC_ZP'] = fdelta_zeropoint

    return  lf
