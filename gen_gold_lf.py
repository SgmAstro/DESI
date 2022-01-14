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
from   delta8_limits import dd8_limits
from   renormalise_d8LF import lumfn_d8_normalise


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

    Area = gama_zmax.meta['AREA']

    if field != None:
        assert  len(found_fields) == 1, 'ERROR: EXPECTED SINGLE FIELD RESTRICTED INPUT, e.g. G9.'

        Area /= 3.
        
    print('Retrieved area {} [sq. deg.]'.format(Area))

    VV = volcom(zmax, Area) - volcom(zmin, Area)
    
    gama_vmax = vmaxer(gama_zmax, zmin, zmax, Area, extra_cols=['MCOLOR_0P0'])

    print('WARNING:  Found {:.3f}% with zmax < 0.0'.format(100. * np.mean(gama_vmax['ZMAX'] <= 0.0)))
    
    # TODO: Why do we need this?                                                                                                   
    gama_vmax = gama_vmax[gama_vmax['ZMAX'] >= 0.0]

    gama_vmax.meta = gama_zmax.meta
    gama_vmax.meta.update({'FORCE_ZMIN': zmin, 'FORCE_ZMAX': zmax, 'VOLUME': VV})
    
    print('Writing {}.'.format(opath))

    gama_vmax.write(opath, format='fits', overwrite=True)
    
    ##  Luminosity fn.
    opath = opath.replace('vmax', 'lumfn')
    
    result = lumfn(gama_vmax, VV)

    result.meta = gama_vmax.meta

    print('Writing {}.'.format(opath))
    
    result.write(opath, format='fits', overwrite=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Generate Gold luminosity function.')
    parser.add_argument('-f', '--field', type=str, help='select equatorial GAMA field: G9, G12, G15', required=False, default=None)
    parser.add_argument('-d', '--density_split', type=bool, help='Trigger density split luminosity function.', default=False)
    parser.add_argument('--dryrun', action='store_true', help='dryrun.')
    
    args = parser.parse_args()

    field = args.field
    dryrun = args.dryrun
    density_split = args.density_split

    print(field, dryrun, density_split)

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
        
        rand_path = os.environ['RANDOMS_DIR'] + '{}/randoms_bd_ddp_n8_{}_0.fits'.format(os.environ['CSCRATCH'], field)
        rand = Table.read(rand_path)
        
        print('Read {}'.format(rand_path))
        
        fpath = os.environ['GOLD_DIR'] + '/gama_gold_zmax.fits'
        
        for idx in range(4):
            ddp_idx   = idx + 1
            ddp_fpath = fpath.replace('zmax', '{}_ddp_n8_d0_{:d}'.format(field, idx))
            ddp_opath = ddp_fpath.split('.')[0] + '_vmax.fits'

            print()
            print('Reading: {}'.format(ddp_fpath))
            
            process_cat(ddp_fpath, ddp_opath, field=field)
        
            print('PROCESS CAT FINISHED.')
                    
            lumfn_path = os.environ['GOLD_DIR'] + '/gama_gold_{}_ddp_n8_d0_{}_lumfn.fits'.format(field, idx)

            print('Reading: {}'.format(lumfn_path))
            
            result = Table.read(lumfn_path)        

            result.pprint()

            scale = rand.meta['DDP1_d{}_VOLFRAC'.format(idx)]

            print('Found d8 renormalisation scale of {:.3f}'.format(scale))
            
            result = lumfn_d8_normalise(result, 1. / scale)

            result.pprint()
            
            lims = dd8_limits[idx]
            d8   = np.mean(lims)

            sc   = named_schechter(result['MEDIAN_M'], named_type='TMR')
            sc  *= (1. + d8) / (1. + 0.007)
            
            result['D8_REFSCH'] = sc 

            print('Writing {}'.format(lumfn_path))
            
            result.write(lumfn_path, format='fits', overwrite=True)

    print('Done.')
