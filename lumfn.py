import numpy as np

from   astropy.table import Table
from   cosmo import volcom

def lumfn(dat, vol, Ms=np.arange(-25.5, -15.5, 0.3), Mcol='MCOLOR_0P0'):
    idxs = np.digitize(dat[Mcol], bins=Ms)

    result = []
    
    for idx in np.unique(idxs):
        sample = dat[idxs == idx]
        nsample = len(sample)
        
        median = np.median(sample[Mcol])
    
        result.append([median, nsample / vol, np.sqrt(nsample) / vol, np.sum(1./sample['VMAX'])])
            
    return  np.array(result)
