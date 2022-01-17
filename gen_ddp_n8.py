import os
import fitsio
import argparse
import numpy as np
import matplotlib.pyplot as plt

from   astropy.table import Table
from   scipy.spatial import KDTree
from   cartesian import cartesian
from   delta8_limits import delta8_tier
from   gama_limits import gama_field


parser = argparse.ArgumentParser(description='Generate DDP1 N8 for all gold galaxies.')
parser.add_argument('-f', '--field', help='GAMA field for randoms_bd file (fill factor & bound_dist retrieval).', required=True)
parser.add_argument('-d', '--dryrun', help='Dryrun.', action='store_true')

args   = parser.parse_args()
field  = args.field.upper()
dryrun = args.dryrun

fpath  = os.environ['GOLD_DIR'] + '/gama_gold_ddp.fits'

if dryrun:
    fpath = fpath.replace('.fits', '_dryrun.fits')

print('Reading: {}'.format(fpath))

# Read ddp cat.    
dat    = Table.read(fpath)

assert 'DDP1_DENS' in dat.meta

points       = np.c_[dat['CARTESIAN_X'], dat['CARTESIAN_Y'], dat['CARTESIAN_Z']]
points       = np.array(points, copy=True)

kd_tree_all  = KDTree(points)

# Read randoms bound_dist.
realz  = 0

rpath  = os.environ['RANDOMS_DIR'] + '/randoms_bd_{}_{:d}.fits'.format(field, realz)

if dryrun:
    rpath = rpath.replace('.fits', '_dryrun.fits')

print('Reading: {}'.format(rpath))

if not os.path.isfile(rpath):
    raise RuntimeError('Expect random bound dist. file for field {}; Run bound_dist.py for this field'.format(field))
    
rand, rand_hdr = fitsio.read(rpath, header=True)

rpoints = np.c_[rand['CARTESIAN_X'], rand['CARTESIAN_Y'], rand['CARTESIAN_Z']]
rpoints = np.array(rpoints, copy=True)

print('Creating big rand. tree.')

big_tree = KDTree(rpoints)

print('Querying tree for closest rand.')

dd, ii   = big_tree.query([x for x in points], k=1)

# Find closest random for bound_dist and fill factor. 
# These randoms are split by field.
dat['RANDSEP']    = dd
dat['RANDMATCH']  = rand['RANDID'][ii]
dat['BOUND_DIST'] = rand['BOUND_DIST'][ii]
dat['FILLFACTOR'] = rand['FILLFACTOR'][ii]

for idx in range(3):
    # Calculate DDP1/2/3 N8 for all gold galaxies.
    ddp_idx      = idx + 1
    
    ddp          = dat[dat['DDP'][:,idx] == 1]
    points_ddp   = np.c_[ddp['CARTESIAN_X'], ddp['CARTESIAN_Y'], ddp['CARTESIAN_Z']]
    points_ddp   = np.array(points_ddp, copy=True)

    print('Building tree for DDP {}'.format(ddp_idx))
    
    kd_tree_ddp  = KDTree(points_ddp)

    print('Querying tree for DDP {}'.format(ddp_idx))

    indexes_ddp  = kd_tree_all.query_ball_tree(kd_tree_ddp, r=8.)

    dat['DDP{:d}_N8'.format(ddp_idx)] = np.array([len(idx) for idx in indexes_ddp])

dat.pprint()

##  Derived.
dat.meta['VOL8']   = (4./3.)*np.pi*(8.**3.)

dat['DDP1_DELTA8'] = (dat['DDP1_N8'] / (dat.meta['VOL8'] * dat.meta['DDP1_DENS']) / dat['FILLFACTOR']) - 1. 
dat['DDP2_DELTA8'] = (dat['DDP2_N8'] / (dat.meta['VOL8'] * dat.meta['DDP2_DENS']) / dat['FILLFACTOR']) - 1. 
dat['DDP3_DELTA8'] = (dat['DDP3_N8'] / (dat.meta['VOL8'] * dat.meta['DDP3_DENS']) / dat['FILLFACTOR']) - 1. 

for x in dat.meta.keys():
    print('{}\t\t{}'.format(x.ljust(20), dat.meta[x]))

print('Writing {}'.format(fpath.replace('ddp', 'ddp_n8')))

dat.write(fpath.replace('ddp', 'ddp_n8'), overwrite=True, format='fits')

# Generate ddp_n8_d0 files for LF(d8) files, limited to DDP1 (and redshift range).
dat = dat[(dat['ZGAMA'] > dat.meta['DDP1_ZMIN']) & (dat['ZGAMA'] < dat.meta['DDP1_ZMAX'])]
dat['DDP1_DELTA8_TIER'] = delta8_tier(dat['DDP1_DELTA8'])

utiers = np.unique(dat['DDP1_DELTA8_TIER'].data)

if -99 in utiers:
    utiers = utiers.tolist()    
    utiers.remove(-99)
    utiers = np.array(utiers)

print(utiers)

if not np.all(utiers == np.arange(4)):
    print('WARNING: MISSING d8 TIERS ({})'.format(utiers))
    
else:
    print(utiers)

print('Delta8 spans {:.4f} to {:.4f} over {} tiers.'.format(dat['DDP1_DELTA8'].min(), dat['DDP1_DELTA8'].max(), utiers))

for tier in utiers:
    print()
    print('---- d{} ----'.format(tier))

    isin  = (dat['DDP1_DELTA8_TIER'].data == tier)    
    to_write = dat[isin]

    assert 'AREA' in dat.meta.keys()
    assert 'AREA' in to_write.meta.keys()
    
    #  E.g. /global/cscratch1/sd/mjwilson/norberg//GAMA4/gama_gold_G9_ddp_n8_d0_0.fits                                                  
    opath = fpath.replace('ddp', 'ddp_n8_d0_{:d}'.format(tier))

    print('Writing {}.'.format(opath))
        
    to_write.write(opath, format='fits', overwrite=True)

    ## 
    isin = to_write['FIELD'] == field
    to_write_field = to_write[isin]
    
    opath_field = opath.replace('gold', 'gold_{}'.format(field))

    print('Writing {}.'.format(opath_field))
    
    to_write_field.write(opath_field, format='fits', overwrite=True)

print('\n\nDone.\n\n')
