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

ngal=1500
Area = 180.
dryrun=False

fpath = os.environ['CSCRATCH'] + '/norberg//GAMA4/gama_gold_zmax.fits'

if dryrun:
    fpath = fpath.replace('_zmax', '_zmax_{:d}k.fits'.format(np.int(ngal / 1000.)))

opath = fpath.replace('zmax', 'vmax')

gama_zmax = Table.read(fpath)

zmin = gama_zmax['ZGAMA'].min()
zmax = gama_zmax['ZGAMA'].max()

gama_vmax = vmaxer(gama_zmax, zmin, zmax, Area, extra_cols=['MCOLOR_0P0'])

# TODO: Why do we need this?
gama_vmax = gama_vmax[gama_vmax['ZMAX'] > 0.0]

gama_vmax.meta = {'FORCE_ZMIN': zmin, 'FORCE_ZMAX': zmax, 'Area': Area}

print('Writing {}.'.format(opath))

gama_vmax.write(opath, format='fits', overwrite=True)


opath = fpath.replace('zmax', 'lumfn')

VV      = volcom(gama_vmax['ZGAMA'].max(), Area) - volcom(gama_vmax['ZGAMA'].min(), Area)
result  = lumfn(gama_vmax, VV)

result.meta = {'FORCE_ZMIN': zmin, 'FORCE_ZMAX': zmax, 'Area': Area, 'Vol': VV}

print('Writing {}.'.format(opath))

result.write(opath, format='fits', overwrite=True)

'''
Ms      = np.arange(-23., -15., 0.01)
tmr_phi = schechter(Ms, named_type='TMR')

pl.plot(Ms, np.log10(tmr_phi), label='TMR')

pl.plot(result[:,0], np.log10(result[:,1]), label='No weights', alpha=0.4)
pl.plot(result[:,0], np.log10(result[:,3]), label='IVMAX')

pl.xlabel(r'$M$')
pl.ylabel(r'$\Phi(M)$')

pl.ylim(-4.25, -1.)
pl.xlim(-23., -15.)

pl.legend(frameon=False, loc=2)
'''
