import numpy as np

from   astropy.table import Table
from   cosmo import volcom

def lumfn(dat, vol, Ms=np.arange(-25.5, -15.5, 0.3), Mcol='MCOLOR_0P0'):
    idxs = np.digitize(dat[Mcol], bins=Ms)

    result = []

    ds = np.diff(Ms)
    dM = ds[0]

    assert np.all(ds == dM)
    
    for idx in np.unique(idxs):
        sample = dat[idxs == idx]
        nsample = len(sample)
        
        median = np.median(sample[Mcol])

        ivmax = 1./sample['VMAX'].data
        variance = np.sum(ivmax**2.) / np.sum(ivmax)**2.
        
        result.append([median, nsample / dM / vol, np.sqrt(nsample) / dM / vol, np.sum(ivmax) / dM, np.sqrt(variance) / dM / vol])

    # TODO: PHI_IVMAX_ERROR is not right.
    return  Table(np.array(result), names=['MEDIAN_M', 'PHI_N', 'PHI_N_ERROR', 'PHI_IVMAX', 'PHI_IVMAX_ERROR'])
