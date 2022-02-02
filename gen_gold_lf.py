import os
import sys
import runtime
import argparse
import pylab as pl
import numpy as np
import astropy.io.fits as fits
import matplotlib.pyplot as plt

from   astropy.table import Table, vstack
from   vmaxer import vmaxer
from   smith_kcorr import test_plots, test_nonnative_plots
from   cosmo import distmod, volcom
from   lumfn import lumfn
from   schechter import schechter, named_schechter
from   gama_limits import gama_field, gama_limits
from   renormalise_d8LF import renormalise_d8LF
from   delta8_limits import d8_limits


def process_cat(fpath, vmax_opath, field=None, rand_paths=[]):
    assert 'vmax' in vmax_opath

    opath = vmax_opath

    if not os.path.isfile(fpath):
        print('WARNING:  Failed to find {}'.format(fpath))
        return 1

    gama_zmax = Table.read(fpath)

    if 'FIELD' not in gama_zmax.dtype.names:
        # print('WARNING:  Missing FIELD keyword, adding it.')
        # gama_zmax['FIELD'] = gama_field(gama_zmax['RA'].data, gama_zmax['DEC'].data)

        raise  RuntimeError('FIELD MISSING FROM DTYPES.')
        
    found_fields = np.unique(gama_zmax['FIELD'].data)
        
    print('Found fields: {}'.format(found_fields))
                            
    zmin = gama_zmax['ZGAMA'].min()
    zmax = gama_zmax['ZGAMA'].max()
    
    print('Found redshift limits: {:.3f} < z < {:.3f}'.format(zmin, zmax))

    if field != None:
        assert  len(found_fields) == 1, 'ERROR: EXPECTED SINGLE FIELD RESTRICTED INPUT, e.g. G9.'

    if len(rand_paths) > 0:
        rand  = vstack([Table.read(_x) for _x in rand_paths])
    else:
        rand  = None

    gama_vmax = vmaxer(gama_zmax, zmin, zmax, extra_cols=['MCOLOR_0P0', 'DDPMALL_0P0_VISZ'], rand=rand)

    print('WARNING:  Found {:.3f}% with zmax < 0.0'.format(100. * np.mean(gama_vmax['ZMAX'] <= 0.0)))
    
    # TODO: Why do we need this?                                                                                                   
    gama_vmax = gama_vmax[gama_vmax['ZMAX'] >= 0.0]
    
    print('Writing {}.'.format(opath))

    gama_vmax.write(opath, format='fits', overwrite=True)
    
    ##  Luminosity fn.
    opath  = opath.replace('vmax', 'lumfn')
    result = lumfn(gama_vmax)

    print('Writing {}.'.format(opath))
    
    result.write(opath, format='fits', overwrite=True)

    return 0


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Generate Gold luminosity function.')
    parser.add_argument('-f', '--field', type=str, help='Select equatorial GAMA field: G9, G12, G15', required=False, default=None)
    parser.add_argument('-d', '--density_split', help='Trigger density split luminosity function.', action='store_true')
    parser.add_argument('--dryrun', action='store_true', help='dryrun.')
    parser.add_argument('--prefix', help='filename prefix', default='randoms')

    parser.add_argument('--nooverwrite',  help='Do not overwrite outputs if on disk', action='store_true')
    
    args   = parser.parse_args()

    field  = args.field
    dryrun = args.dryrun
    density_split = args.density_split
    prefix = args.prefix

    if density_split:
        assert  field != None
        assert  'ddp1' in prefix

    if not density_split:
        print('Generating Gold reference LF.')

        # Bounded by gama gold, reference schechter limits:  
        # 0.039 < z < 0.263.
        # Note: not split by field. 
        fpath = os.environ['GOLD_DIR'] + '/gama_gold_ddp.fits'
        
        if dryrun:
            fpath = fpath.replace('.fits', '_dryrun.fits')

        opath = fpath.replace('ddp', 'vmax')

        # 
        if args.nooverwrite:
            if os.path.isfile(opath) and os.path.isfile(opath.replace('vmax', 'lumfn')):
                
                print('{} found on disk and overwrite forbidden (--nooverwrite).'.format(fpath))
                exit(0)
            
        '''
        all_rpaths = [os.environ['RANDOMS_DIR'] + '/{}_bd_ddp_n8_G{}_0.fits'.format(prefix, ff) for ff in [9, 12, 15]]

        if dryrun:
            all_rpaths = [_rpath.replace('.fits', '_dryrun.fits') for _rpath in all_rpaths]
        '''

        print(f'Reading: {fpath}')
        print(f'Writing: {opath}')

        process_cat(fpath, opath, rand_paths=[])

    else:
        print('Generating Gold density-split LF.')

        field = field.upper()

        rpath = os.environ['RANDOMS_DIR'] + '/{}_bd_ddp_n8_{}_0.fits'.format(prefix, field)
        
        if dryrun:
            rpath = rpath.replace('.fits', '_dryrun.fits')
                
        if dryrun:
            # A few galaxies have a high probability to be in highest density only. 
            utiers = np.array([8])

        else:
            utiers = np.arange(len(d8_limits))
                    
        all_rands = None 

        for idx in utiers:
            ddp_idx   = idx + 1

            # Bounded by DDP1 z limits. 
            ddp_fpath = os.environ['GOLD_DIR'] + '/gama_gold_{}_ddp_n8_d0_{:d}.fits'.format(field, idx)
            ddp_opath = ddp_fpath.split('.')[0] + '_vmax.fits'
            
            if dryrun:
                ddp_fpath = ddp_fpath.replace('.fits', '_dryrun.fits')
                ddp_opath = ddp_opath.replace('.fits', '_dryrun.fits')

            print()
            print('Reading: {}'.format(ddp_fpath))
            
            failure = process_cat(ddp_fpath, ddp_opath, field=field, rand_paths=[rpath])

            if failure:
                print('FAILED on d0 tier {:d}; skipping.'.format(idx))

                continue
        
            print('LF process cat. complete.')
                    
            result = Table.read(ddp_opath.replace('vmax', 'lumfn'))        

            # result.pprint()

            if all_rands == None:
                all_rpaths = [os.environ['RANDOMS_DIR'] + '/{}_bd_ddp_n8_G{}_0.fits'.format(prefix, ff) for ff in [9, 12, 15]]

                if dryrun:
                    all_rpaths = [_rpath.replace('.fits', '_dryrun.fits') for _rpath in all_rpaths]

                all_rands = [Table.read(_x) for _x in all_rpaths]

            # Calculated for DDP1 redshift limits. 
            fdelta = np.array([x.meta['DDP1_d{}_VOLFRAC'.format(idx)] for x in all_rands])
            d8     = np.array([x.meta['DDP1_d{}_TIERMEDd8'.format(idx)] for x in all_rands])

            print('Field vol renormalization: {}'.format(fdelta))
            print('Field d8  renormalization: {}'.format(d8))

            fdelta = fdelta.mean()
            d8     = d8.mean()

            print('Found mean vol. renormalisation scale of {:.3f}'.format(fdelta))
            print('Found mean  d8  renormalisation scale of {:.3f}'.format(d8))
            
            result = renormalise_d8LF(result, fdelta)
            
            result['REF_SCHECHTER']  = named_schechter(result['MEDIAN_M'], named_type='TMR')
            result['REF_SCHECHTER'] *= (1. + d8) / (1. + 0.007)

            print('LF renormalization and ref. schechter complete.')
            
            result.pprint()

            # 
            sch_Ms = np.arange(-23., -15., 1.e-3)

            sch    = named_schechter(sch_Ms, named_type='TMR')
            sch   *= (1. + d8) / (1. + 0.007)

            ##
            ref_result = Table(np.c_[sch_Ms, sch], names=['MS', 'd{}_REFSCHECHTER'.format(idx)])            
            ref_result.meta['DDP1_d{}_VOLFRAC'.format(idx)]   = fdelta
            ref_result.meta['DDP1_d{}_TIERMEDd8'.format(idx)] = d8

            print('Writing {}'.format(ddp_opath.replace('vmax', 'lumfn')))

            primary_hdu    = fits.PrimaryHDU()

            keys           = sorted(result.meta.keys())
            
            header         = {}
            
            for key in keys:
                header[key] = str(result.meta[key])

            hdr            = fits.Header(header)
            result_hdu     = fits.BinTableHDU(result, name='LUMFN', header=hdr)
            ref_result_hdu = fits.BinTableHDU(ref_result, name='REFERENCE')
            hdul           = fits.HDUList([primary_hdu, result_hdu, ref_result_hdu])

            hdul.writeto(ddp_opath.replace('vmax', 'lumfn'), overwrite=True, checksum=True)

    print('Done.')
