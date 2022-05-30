import fitsio
import subprocess
import astropy.io.fits as      fits
import numpy           as      np

from   astropy.table   import  Table
from   cosmo           import  volcom
from   schechter       import  named_schechter


def multifield_lumfn(lumfn_list, ext=None, weight=None, sub_cols=None):
    if ext is None:
        tables = [Table.read(x) for x in lumfn_list]
    else:
        tables = [Table.read(x, ext) for x in lumfn_list]

    if weight is not None:
        weights = np.array([tab.meta[weight] for tab in tables]).astype(float)

        print('Retrieved relative weights: {} for {} weight.'.format(weights / np.sum(weights), weight))

    else:
        weights = None

    def sum_rule(tables, col, weights=None):
        data = [table[col].data for table in tables]
        data = np.c_[data].T
        
        return np.sum(data, axis=1)

    def mean_rule(tables, col, weights=None):
        data = [table[col].data for table in tables]
        data = np.c_[data].T
        
        return np.average(data, axis=1, weights=weights)

    def quadsum_rule(tables, col, weights=None):
        data = [table[col].data for table in tables]
        data = np.c_[data].T
        
        # TODO
        if weights is not None:
            pass 

            # print('WARNING: weights is unsupported for lumfn quadsum rule.')

        return  np.sqrt(np.sum(data**2., axis=1))

    
    result = Table()

    if ext in [None, 'LUMFN']:
        sum_cols   = ['N']
        mean_cols  = ['MEDIAN_M', 'MEAN_M', 'MID_M', 'IVMAXMEAN_M', 'PHI_N', 'PHI_IVMAX', 'V_ON_VMAX', 'REF_SCHECHTER', 'REF_RATIO', 'PHI_STEPWISE']
        qsum_cols  = ['PHI_N_ERROR', 'PHI_IVMAX_ERROR']

    elif ext == 'LUMFN_STEP':
        sum_cols   = ['N']
        mean_cols  = ['MID_M', 'PHI_STEPWISE', 'REF_RATIO']
        qsum_cols  = []
        
    elif ext == 'REFERENCE':
        sum_cols   = []
        mean_cols  = ['MS', 'REFSCHECHTER']
        qsum_cols  = [] 

    else:
        raise  RuntimeError(f'MultifieldLumfn:  Extension {ext} is not supported.')

    if sub_cols != None:
        sum_cols  = [x for x in sum_cols  if x in sub_cols]
        mean_cols = [x for x in mean_cols if x in sub_cols]
        qsum_cols = [x for x in qsum_cols if x in sub_cols]

    for m in mean_cols:
        result[m] = mean_rule(tables, m, weights=weights)

    for s in sum_cols:
        result[s] = sum_rule(tables, s, weights=weights)
        
    for q in qsum_cols:
        result[q] = quadsum_rule(tables, q, weights=weights)

    if ext != 'REFERENCE':
        result['VALID'] = result['N'] >= 5
    
    return  result

def lumfn(dat, Ms=None, Mcol='MCOLOR_0P0', jackknife=None, opath=None, d8=None):
    if type(jackknife) == np.ndarray:
        for jk in jackknife:
            lumfn(dat, Ms=Ms, Mcol=Mcol, jackknife=int(jk), opath=opath)

        return 0
    
    elif type(jackknife) == int:
        pass

    elif jackknife is None:
        pass

    else:
        raise ValueError('Unsupported jackknife of type {}'.format(type(jackknife)))

    if Ms == None:
        # np.arange(-25.5, -15.5, 0.2) 
        Ms = np.linspace(-23.,  -16.,  36)
                
    dat   = Table(dat, copy=True)

    # If values in x are beyond the bounds of bins, 0 or len(bins) is returned as appropriate.
    keep  = (dat[Mcol] >= Ms.min()) & (dat[Mcol] <= Ms.max())
    dat   = dat[keep]

    dvmax = dat['VMAX'].data
    vol   = dat.meta['FORCE_VOL']
    
    # default:  bins[i-1] <= x < bins[i]
    
    if jackknife is not None:
        print('Solving for jack knife {}'.format(jackknife))

        jk_volfrac = dat.meta['JK_VOLFRAC']

        vol       *= jk_volfrac

        dat        = dat[dat['JK'] != f'JK{jackknife}']
        dvmax      = jk_volfrac * dat['VMAX'].data
    
    idxs   = np.digitize(dat[Mcol], bins=Ms)
    result = []

    print('\n\nSolving for Ms: {}'.format(Ms))

    ds     = np.diff(Ms)

    ds     = np.round(ds, decimals=4)
    dM     = ds[0]

    assert  np.all(ds == dM)
    
    for ii, idx in enumerate(np.arange(0, len(Ms), 1)):
        sample  = dat[idxs == idx]
        nsample = len(sample)

        # print(sample)
        
        vmax    = dvmax[idxs == idx]

        ivmax   = 1. / vmax
        ivmax2  = 1. / vmax**2.

        if nsample > 0:
            median = np.median(sample[Mcol])
            mean   = np.mean(sample[Mcol])
            wmean  = np.average(sample[Mcol], weights=ivmax) 
            mid    = Ms[ii] + dM/2.

        else:
            median = Ms[ii] + dM/2.
            mean   = median
            wmean  = mean
            mid    = mean

        # print(median)

        if len(vmax) == 0:
            median_vmax = 0

        else:
            median_vmax = np.median(vmax) / vol

        result.append([median,\
                       mean,\
                       mid,\
                       wmean,\
                       nsample / dM / vol,\
                       np.sqrt(nsample) / dM / vol,\
                       np.sum(ivmax) / dM,\
                       np.sqrt(np.sum(ivmax2)) / dM,\
                       nsample,
                       median_vmax])

    names  = ['MEDIAN_M', 'MEAN_M', 'MID_M', 'IVMAXMEAN_M', 'PHI_N', 'PHI_N_ERROR', 'PHI_IVMAX', 'PHI_IVMAX_ERROR', 'N', 'V_ON_VMAX']

    result = Table(np.array(result), names=names)
    result['VALID'] = result['N'] >= 5.
    result['REF_SCHECHTER']       = named_schechter(result['MEDIAN_M'], named_type='TMR')

    if d8 != None:
        # TODO HARDCODE 0.007                                                                                                                                                                             
        result['REF_SCHECHTER']  *= (1. + d8) / (1. + 0.007)
        result.meta['DDP1_D8']    = d8     

    result['REF_RATIO']           = result['PHI_IVMAX'] / result['REF_SCHECHTER']

    result.meta.update(dat.meta)

    result.pprint()
    
    result.meta['MS']             = str(['{:.4f}'.format(x) for x in Ms.tolist()])
    result.meta['FORCE_VOL']      = vol
    result.meta['ABSMAG_DEF']     = Mcol
    result.meta['EXTNAME']        = 'LUMFN'
    
    if jackknife is not None:        
        result.meta['EXTNAME']    = 'LUMFN_JK{}'.format(jackknife)
        result.meta['RENORM']     = 'FALSE'
        result.meta['JK_VOLFRAC'] = dat.meta['JK_VOLFRAC']
        result.meta['NJACK']      = dat.meta['NJACK']
        result                    = fits.convenience.table_to_hdu(result)

        with fits.open(opath, mode='update') as hdulist:
            hdulist.append(result)
            hdulist.flush()  
            hdulist.close()

        cmds   = []
        cmds.append(f'chgrp desi {opath}')
        cmds.append(f'chmod  700 {opath}')

        for cmd in cmds:
            output = subprocess.check_output(cmd, shell=True)

            print(cmd, output)

        return  0

    else:
        return  result 
