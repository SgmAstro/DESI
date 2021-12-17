import os
import fitsio
from astropy.table import Table
from ddp import get_ddps


fpath = os.environ['CSCRATCH'] + '/norberg/GAMA4/gama_gold_zmax.fits'
dat   = fitsio.read(fpath)
dat   = Table(dat)

Area  = 180. 

dat['DDP'], zlims = get_ddps(Area, dat['MCOLOR_0P0'], dat['ZGAMA'])
dat.meta = zlims

print(zlims)

dat.write(fpath.replace('zmax', 'ddp'), format='fits', overwrite=True)
