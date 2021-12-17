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
from   schechter import schechter
from   gama_limits import gama_field

def process_cat(fpath, vmax_opath):
    assert 'vmax' in vmax_opath

    opath = vmax_opath

    gama_zmax = Table.read(fpath)

    # HACK:  remove after pipeline rerun with G9, G12, G15 randoms.
    print('HACK: Limiting to G9.')
        
    gama_zmax['FIELD'] = gama_field(gama_zmax['RA'], gama_zmax['DEC'])
    gama_zmax = gama_zmax[gama_zmax['FIELD'] == 'G9']
        
    zmin = gama_zmax['ZGAMA'].min()
    zmax = gama_zmax['ZGAMA'].max()

    gama_vmax = vmaxer(gama_zmax, zmin, zmax, Area, extra_cols=['MCOLOR_0P0'])

    # TODO: Why do we need this?                                                                                                                                    
    gama_vmax = gama_vmax[gama_vmax['ZMAX'] > 0.0]

    gama_vmax.meta = {'FORCE_ZMIN': zmin, 'FORCE_ZMAX': zmax, 'Area': Area}

    print('Writing {}.'.format(opath))

    gama_vmax.write(opath, format='fits', overwrite=True)

    ##  Luminosity fn.
    opath = vmax_opath.replace('vmax', 'lumfn')

    VV      = volcom(gama_vmax['ZGAMA'].max(), Area) - volcom(gama_vmax['ZGAMA'].min(), Area)
    result  = lumfn(gama_vmax, VV)

    result.meta = {'FORCE_ZMIN': zmin, 'FORCE_ZMAX': zmax, 'Area': Area, 'Vol': VV}

    print('Writing {}.'.format(opath))

    result.write(opath, format='fits', overwrite=True)
    
    
ngal=1500
Area = 180.
dryrun=False
density_split=True

if not density_split:
    fpath = os.environ['CSCRATCH'] + '/norberg//GAMA4/gama_gold_zmax.fits'

    if dryrun:
        fpath = fpath.replace('_zmax', '_zmax_{:d}k.fits'.format(np.int(ngal / 1000.)))

    opath = fpath.replace('zmax', 'vmax')
        
    process_cat(fpath, opath)

else:
    fpath = os.environ['CSCRATCH'] + '/norberg//GAMA4/gama_gold_zmax.fits'

    for idx in range(4):
        ddp_idx   = idx + 1
        ddp_fpath = fpath.replace('zmax', 'ddp_n8_d0_{:d}'.format(idx))
        ddp_opath = ddp_fpath.split('.')[0] + '_vmax.fits'

        print()
        print(ddp_fpath)
        print(ddp_opath)
        
        process_cat(ddp_fpath, ddp_opath)
