import os
import argparse
import fitsio

from   astropy.table import Table
from   ddp import get_ddps


parser = argparse.ArgumentParser(description='Gen ddp cat.')
parser.add_argument('-d', '--dryrun', help='Dryrun.', action='store_true')

args = parser.parse_args()
dryrun = args.dryrun

fpath = os.environ['GOLD_DIR'] + '/gama_gold_zmax.fits'

if dryrun:
    fpath = fpath.replace('.fits', '_dryrun.fits')

print('Reading: {}'.format(fpath))
    
dat   = fitsio.read(fpath)
dat   = Table(dat)

Area  = dat.meta['AREA']

print('Judging DDP.')

dat['DDP'], zlims = get_ddps(Area, dat['MCOLOR_0P0'], dat['ZGAMA'])
dat.meta.update(zlims)

print(Area)
print(zlims)

opath = fpath.replace('zmax', 'ddp')

print('Reading: {}'.format(opath))

dat.write(opath, format='fits', overwrite=True)
