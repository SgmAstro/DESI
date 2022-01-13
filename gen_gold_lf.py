import os
import sys
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

import argparse

def process_cat(fpath, vmax_opath, field=None):
    assert 'vmax' in vmax_opath

    opath = vmax_opath

    gama_zmax = Table.read(fpath)

    if 'FIELD' not in gama_zmax.dtype.names:
        gama_zmax['FIELD'] = gama_field(gama_zmax['RA'].data, gama_zmax['DEC'].data)
    
    if field != None:
        assert field in gama_limits.keys()
                
        gama_zmax = gama_zmax[gama_zmax['FIELD'].data == field]
                    
    zmin = gama_zmax['ZGAMA'].min()
    zmax = gama_zmax['ZGAMA'].max()

    gama_vmax = vmaxer(gama_zmax, zmin, zmax, Area, extra_cols=['MCOLOR_0P0'])

    # TODO: Why do we need this?                                                                                                                                    
    gama_vmax = gama_vmax[gama_vmax['ZMAX'] > 0.0]

    gama_vmax.meta = {'FORCE_ZMIN': zmin, 'FORCE_ZMAX': zmax, 'Area': Area}

    print('Writing {}.'.format(opath))

    gama_vmax.write(opath, format='fits', overwrite=True)

    ##  Luminosity fn.
    opath = opath.replace('vmax', 'lumfn')

    VV = volcom(gama_vmax['ZGAMA'].max(), Area) - volcom(gama_vmax['ZGAMA'].min(), Area)

    result = lumfn(gama_vmax, VV)

    result.meta = {'FORCE_ZMIN': zmin, 'FORCE_ZMAX': zmax, 'Area': Area, 'Vol': VV}

    print('Writing {}.'.format(opath))
    
    result.write(opath, format='fits', overwrite=True)
    
ngal = 1500
Area = 180.
dryrun=False
density_split=True

parser = argparse.ArgumentParser(description='Select GAMA field.')
parser.add_argument('-f', '--field', type=str, help='select equatorial GAMA field: G9, G12, G15', required=True)
args = parser.parse_args()
field = args.field.upper()

if not density_split:
    field = ''
    print('IGNORING FIELD GENERATING G9 to G15')
    
    fpath = os.environ['CSCRATCH'] + '/norberg/GAMA4/gama_gold_zmax.fits'

    if dryrun:
        fpath = fpath.replace('_zmax', '_zmax_{:d}k.fits'.format(np.int(ngal / 1000.)))

    opath = fpath.replace('zmax', 'vmax')
        
    process_cat(fpath, opath)

else:
    fpath = os.environ['CSCRATCH'] + '/norberg/GAMA4/gama_gold_zmax.fits'

    for idx in range(4):
        ddp_idx   = idx + 1
        ddp_fpath = fpath.replace('zmax', '{}_ddp_n8_d0_{:d}'.format(field, idx))
        ddp_opath = ddp_fpath.split('.')[0] + '_vmax.fits'

        print()
        print(ddp_fpath)
        print(ddp_opath)
        
        process_cat(ddp_fpath, ddp_opath, field=field)
        
        print('PROCESS CAT FINISHED.')
        
        rand_path = '{}/desi/BGS/Sam/randoms_bd_ddp_n8_{}_0.fits'.format(os.environ['CSCRATCH'], field)
        rand = Table.read(rand_path)
        scale = rand.meta['DDP1_d{}_VOLFRAC'.format(idx)]
        
        lumfn_path = os.environ['CSCRATCH'] + '/norberg/GAMA4/gama_gold_{}_ddp_n8_d0_{}_lumfn.fits'.format(field, idx)
        result = Table.read(lumfn_path)        
        result = lumfn_d8_normalise(result, scale)
                
        gama_lf_path = os.environ['CSCRATCH'] + '/norberg/GAMA4/gama_gold_lumfn.fits'
        gama_lf =  Table.read(gama_lf_path)
        sc = named_schechter(gama_lf['MEDIAN_M'], named_type='TMR')
        lims = dd8_limits[idx]
        d8 = np.mean(lims)
        sc *= (1. + d8) / (1. + 0.007)
        gama_lf['sc'] = sc 
        
        gama_lf.write(gama_lf_path, format='fits', overwrite=True)

