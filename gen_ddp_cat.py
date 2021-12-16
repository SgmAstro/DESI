import os
import fitsio
from astropy.table import Table
from ddp import ddp

fpath = os.environ['CSCRATCH'] + '/norberg/GAMA4/gama_gold_kE_20k.fits'
dat = fitsio.read(fpath)
dat = Table(dat)

ddp_array = ddp(dat['MCOLOR_0P0'], dat['ZGAMA'])

dat['DDP1'] = ddp_array[0]
dat['DDP2'] = ddp_array[1]
dat['DDP3'] = ddp_array[2]

dat.write(fpath.replace('20k', '20k_ddp'), format='fits', overwrite=True)