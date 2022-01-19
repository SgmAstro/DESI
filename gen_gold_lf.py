import os
import sys
import argparse
import pylab as pl
import numpy as np
import matplotlib.pyplot as plt

from   astropy.table import Table
from   vmaxer import vmaxer
from   smith_kcorr import test_plots, test_nonnative_plots
from   cosmo import distmod, volcom
from   lumfn import lumfn
from   schechter import schechter, named_schechter
from   gama_limits import gama_field, gama_limits
from   renormalise_d8LF import renormalise_d8LF


def process_cat(fpath, vmax_opath, field=None):
    assert 'vmax' in vmax_opath

    opath = vmax_opath

    gama_zmax = Table.read(fpath)

    if 'FIELD' not in gama_zmax.dtype.names:
        print('WARNING:  Missing FIELD keyword, adding it.')
        
        gama_zmax['FIELD'] = gama_field(gama_zmax['RA'].data, gama_zmax['DEC'].data)
        
    found_fields = np.unique(gama_zmax['FIELD'].data)
        
    print('Found fields: {}'.format(found_fields))
                            
    zmin = gama_zmax['ZGAMA'].min()
    zmax = gama_zmax['ZGAMA'].max()
    
    print('Found redshift limits: {:.3f} < z < {:.3f}'.format(zmin, zmax))

    if field != None:
        assert  len(found_fields) == 1, 'ERROR: EXPECTED SINGLE FIELD RESTRICTED INPUT, e.g. G9.'
    
    gama_vmax = vmaxer(gama_zmax, zmin, zmax, extra_cols=['MCOLOR_0P0'])

    print('WARNING:  Found {:.3f}% with zmax < 0.0'.format(100. * np.mean(gama_vmax['ZMAX'] <= 0.0)))
    
    # TODO: Why do we need this?                                                                                                   
    gama_vmax = gama_vmax[gama_vmax['ZMAX'] >= 0.0]
    
    print('Writing {}.'.format(opath))

    gama_vmax.write(opath, format='fits', overwrite=True)
    
    ##  Luminosity fn.
    opath = opath.replace('vmax', 'lumfn')
    
    result = lumfn(gama_vmax)

    print('Writing {}.'.format(opath))
    
    result.write(opath, format='fits', overwrite=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Generate Gold luminosity function.')
    parser.add_argument('-f', '--field', type=str, help='select equatorial GAMA field: G9, G12, G15', required=False, default=None)
    parser.add_argument('-d', '--density_split', help='Trigger density split luminosity function.', action='store_true')
    parser.add_argument('--dryrun', action='store_true', help='dryrun.')
    
    args = parser.parse_args()

    field = args.field
    dryrun = args.dryrun
    density_split = args.density_split
    
    print(field, dryrun, density_split)

    user = os.environ['USER']
    
    if not density_split:
        print('Generating Gold reference LF.')
        
        field = ''
        
        print('IGNORING FIELD ARG., GENERATING ALL OF G9-G15')
    
        fpath = os.environ['GOLD_DIR'] + '/gama_gold_zmax.fits'
        
        if dryrun:
            fpath = fpath.replace('.fits', '_dryrun.fits')

        opath = fpath.replace('zmax', 'vmax')
        
        process_cat(fpath, opath)

    else:
        print('Generating Gold density-split LF.')

        assert field != None

        field = field.upper()
        
        rpath = os.environ['RANDOMS_DIR'] + '/randoms_bd_ddp_n8_{}_0.fits'.format(field)
        
        if dryrun:
            rpath = rpath.replace('.fits', '_dryrun.fits')

        rand = Table.read(rpath)
        
        print('Read {}'.format(rpath))
        
        if dryrun:
            # A few galaxies have a high probability to be in highest density only. 
            utiers = np.array([3])

        else:
            utiers = np.arange(4)
                    
        for idx in utiers:
            ddp_idx   = idx + 1
            ddp_fpath = os.environ['GOLD_DIR'] + '/gama_gold_{}_ddp_n8_d0_{:d}.fits'.format(field, idx)
            ddp_opath = ddp_fpath.split('.')[0] + '_vmax.fits'
            
            if dryrun:
                ddp_fpath = ddp_fpath.replace('.fits', '_dryrun.fits')
                ddp_opath = ddp_opath.replace('.fits', '_dryrun.fits')

            print()
            print('Reading: {}'.format(ddp_fpath))
            
            process_cat(ddp_fpath, ddp_opath, field=field)
        
            print('PROCESS CAT FINISHED.')
                    
            result = Table.read(ddp_opath.replace('vmax', 'lumfn'))        

            result.pprint()
            
            rands = [Table.read(os.environ['RANDOMS_DIR'] + '/randoms_bd_ddp_n8_G{}_0.fits'.format(field)) for field in [9, 12, 15]]

            scale = np.array([x.meta['DDP1_d{}_VOLFRAC'.format(idx)] for x in rands])
            d8    = np.array([x.meta['DDP1_d{}_TIERMEDd8'.format(idx)] for x in rands])

            print('Field vol renormalization: {}'.format(scale))
            print('Field d8  renormalization: {}'.format(d8))

            scale = scale.mean()
            d8    = d8.mean()

            print('Found mean vol. renormalisation scale of {:.3f}'.format(scale))
            print('Found mean  d8  renormalisation scale of {:.3f}'.format(d8))
            
            result = renormalise_d8LF(result, 1. / scale)

            result.pprint()
            
            sc   = named_schechter(result['MEDIAN_M'], named_type='TMR')
            sc  *= (1. + d8) / (1. + 0.007)
            
            result['d{}_REFSCHECHTER'.format(idx)] = sc 

            result.meta['DDP1_d{}_VOLFRAC'.format(idx)]   = scale
            result.meta['DDP1_d{}_TIERMEDd8'.format(idx)] = d8

            print('Writing {}'.format(ddp_opath.replace('vmax', 'lumfn')))
            
            result.write(ddp_opath.replace('vmax', 'lumfn'), format='fits', overwrite=True)

    print('Done.')
