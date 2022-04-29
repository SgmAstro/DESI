import numpy         as     np

from   astropy.table import Table
from   cosmo         import volcom


def multifield_lumfn(lumfn_list, ext=None, weight=None):
    if ext is None:
        tables = [Table.read(x) for x in lumfn_list]
    else:
        tables = [Table.read(x, ext) for x in lumfn_list]

    if weight is not None:
        weights = np.array([tab.meta[weight] for tab in tables]).astype(float)

        print('Retrieved relative weights: {} for {} weight.'.format(weights, weight))

    else:
        weights = None

    def sum_rule(tables, col, weights=None):
        data = [table[col].data for table in tables]
        data = np.c_[data].T
        
        return np.sum(data, axis=1)

    def mean_rule(tables, col, weights=None):
        data = [table[col].data for table in tables]
        data = np.c_[data].T

        if weights is not None:
            print('Mean rule:  applying relative weights.')
        
        return np.average(data, axis=1, weights=weights)

    def quadsum_rule(tables, col, weights=None):
        data = [table[col].data for table in tables]
        data = np.c_[data].T
        
        # TODO 
        if weights is not None:
            print('WARNING: weights is unsupported for lumfn quadsum rule.')

        return  np.sqrt(np.sum(data**2., axis=1))

    
    result    = Table()

    if ext in [None, 'LUMFN']:
        sum_cols   = ['N']
        mean_cols  = ['MEDIAN_M', 'PHI_N', 'PHI_IVMAX', 'V_ON_VMAX', 'REF_SCHECHTER', 'REF_RATIO']
        qsum_cols  = ['PHI_N_ERROR', 'PHI_IVMAX_ERROR']
        
    elif ext == 'REFERENCE':
        sum_cols   = []
        mean_cols  = ['MS', 'REFSCHECHTER']
        qsum_cols  = [] 

    else:
        raise  RuntimeError(f'MultifieldLumfn:  Extension {ext} is not supported.')

    for m in mean_cols:
        result[m] = mean_rule(tables, m, weights=weights)

    for s in sum_cols:
        result[s] = sum_rule(tables, s, weights=weights)
        
    for q in qsum_cols:
        result[q] = quadsum_rule(tables, q, weights=weights)
    
    return  result

def lumfn(dat, Ms=np.arange(-25.5, -15.5, 0.4), Mcol='MCOLOR_0P0', bitmask='IN_D8LUMFN', jk=False, writeto=None):
    dat   = Table(dat, copy=True)

    dat   = dat[dat[bitmask] == 0]

    dvmax = dat['VMAX'].data
    vol   = dat.meta['VOLUME']
    
    # assert  dat[Mcol].min() >= Ms.min()
    # assert  dat[Mcol].max() <= Ms.max()

    # default:  bins[i-1] <= x < bins[i]
    
    if jk:
        n_jk = len(np.unique(dat['JK']))

        dat = dat[dat['JK'] != jk]
        dvmax = dat['VMAX'].data

        vol = vol * (n_jk -1) / n_jk

        # TESTING
        print(n_jk)
        print(vol)
    
    
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
            median = 0.5 * (Ms[idx] + Ms[idx+1])

        vmax    = dvmax[idxs == idx]

        ivmax   = 1. / vmax
        ivmax2  = 1. / vmax**2.

        if len(vmax) == 0:
            median_vmax = 0
        else:
            median_vmax = np.median(vmax) / vol

        result.append([median,\
                       nsample / dM / vol,\
                       np.sqrt(nsample) / dM / vol,\
                       np.sum(ivmax) / dM,\
                       np.sqrt(np.sum(ivmax2)) / dM,\
                       nsample,
                       median_vmax])

        print(result[-1])

    names  = ['MEDIAN_M', 'PHI_N', 'PHI_N_ERROR', 'PHI_IVMAX', 'PHI_IVMAX_ERROR', 'N', 'V_ON_VMAX']

    result = Table(np.array(result), names=names)
    result.meta.update(dat.meta)
    
    result.meta['MS']         = str(['{:.4f}'.format(x) for x in Ms.tolist()])
    result.meta['VOLUME']     = vol
    result.meta['ABSMAG_DEF'] = Mcol

    
    if writeto != None:
        
        opath = writeto

        if jk:
            jk_array = np.unique(dat['JK'])    
            
            for idx in jk_array:
                ##  
                keys           = sorted(result.meta.keys())
                header         = {}

                for key in keys:
                    header[key] = str(result.meta[key])

                primary_hdu    = fits.PrimaryHDU()
                hdr            = fits.Header(header)
                result_hdu     = fits.BinTableHDU(result, name='LUMFN', header=hdr)
                hdul           = fits.HDUList([primary_hdu, result_hdu])

                hdul.writeto(opath, overwrite=True, checksum=True)

            print('Writing {}.'.format(opath))
    
        else:
            result.writeto(opath, overwrite=True, checksum=True)
            
    else:
        print('WARNING: writeto == None')
            
    return  result 
