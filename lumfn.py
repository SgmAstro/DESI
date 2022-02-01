import numpy as np

from   astropy.table import Table
from   cosmo import volcom


def lumfn(dat, Ms=np.arange(-25.5, -15.5, 0.2), Mcol='MCOLOR_0P0'):
    dat = Table(dat, copy=True)

    if 'IN_SAMPLE' in dat.dtype.names:
        # e.g. FILLFACTOR > 0.8 cut   
        dat      = dat[dat['IN_SAMPLE'] > 0]
        vol_frac = dat.meta['IN_SAMPLE_VOLFRAC']

    else:
        vol_frac = 1.

    vol    = dat.meta['VOLUME']
    vol   *= vol_frac

    idxs   = np.digitize(dat[Mcol], bins=Ms)
    result = []

    ds     = np.diff(Ms)
    dM     = ds[0]

    assert  np.all(ds == dM)
    
    for idx in np.arange(len(Ms) - 1):
        sample  = dat[idxs == idx]
        nsample = len(sample)
        
        median  = np.median(sample[Mcol])

        vmax    = sample['VMAX'].data
        vmax   *= vol_frac

        ivmax   = 1. / vmax
        ivmax2  = 1. / vmax**2.

        # TODO: remove NaNs from dataset by setting M to mid bin.
        
        # nsample == 0; set M to mid bin.

        result.append([median,\
                       nsample / dM / vol,\
                       np.sqrt(nsample) / dM / vol,\
                       np.sum(ivmax) / dM,\
                       np.sqrt(np.sum(ivmax2)) / dM,\
                       nsample,
                       np.median(vmax) / vol])

    names = ['MEDIAN_M', 'PHI_N', 'PHI_N_ERROR', 'PHI_IVMAX', 'PHI_IVMAX_ERROR', 'N', 'V_ON_VMAX']

    result = Table(np.array(result), names=names)
    result.meta.update(dat.meta)

    result.meta['MS']     = str(Ms.tolist())
    result.meta['VOLUME'] = vol

    return  result 
