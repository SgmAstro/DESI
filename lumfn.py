import numpy         as     np

from   astropy.table import Table
from   cosmo         import volcom


def multifield_lumfn(lumfn_list):
    tables = [Table.read(x) for x in lumfn_list]
    
    def sum_rule(tables, col):
        data = [table[col].data for table in tables]
        data = np.c_[data].T
        
        return np.sum(data, axis=1)

    def mean_rule(tables, col):
        data = [table[col].data for table in tables]
        data = np.c_[data].T
        
        return np.mean(data, axis=1)

    def quadsum_rule(tables, col):
        data = [table[col].data for table in tables]
        data = np.c_[data].T
        
        return  np.sqrt(np.sum(data**2., axis=1))

    
    result   = Table()
    
    sum_cols  = ['N']
    mean_cols = ['MEDIAN_M', 'PHI_N', 'PHI_IVMAX', 'V_ON_VMAX', 'REF_SCHECHTER']
    qsum_cols = ['PHI_N_ERROR', 'PHI_IVMAX_ERROR']
        
    for m in mean_cols:
        result[m] = mean_rule(tables, m)

    for s in sum_cols:
        result[s] = sum_rule(tables, s)
        
    for q in qsum_cols:
        result[q] = quadsum_rule(tables, q)
    
    return  result

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

    # assert  dat[Mcol].min() >= Ms.min()
    # assert  dat[Mcol].max() <= Ms.max()

    # default:  bins[i-1] <= x < bins[i]
    idxs   = np.digitize(dat[Mcol], bins=Ms)

    result = []

    ds     = np.diff(Ms)
    dM     = ds[0]

    assert  np.all(ds == dM)
    
    for idx in np.arange(len(Ms) - 1):
        sample  = dat[idxs == idx]
        nsample = len(sample)

        if nsample > 0:
            median = np.median(sample[Mcol])
        else:
            # TODO:
            median = 0.5 * (Ms[idx] + Ms[idx+1])

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
    
    result.meta['MS']         = str(Ms.tolist())
    result.meta['VOLUME']     = vol
    result.meta['ABSMAG_DEF'] = Mcol

    return  result 
