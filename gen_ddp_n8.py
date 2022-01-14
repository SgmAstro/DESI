import os
import fitsio
import numpy as np
import matplotlib.pyplot as plt

from   astropy.table import Table
from   scipy.spatial import KDTree
from   cartesian import cartesian
from   delta8_limits import delta8_tier
from   gama_limits import gama_field

import argparse


parser = argparse.ArgumentParser(description='Select GAMA field.')
parser.add_argument('-f', '--field', type=str, help='select equatorial GAMA field: G9, G12, G15', required=True)
args = parser.parse_args()
field = args.field.upper()

fpath = os.environ['CSCRATCH'] + '/norberg/GAMA4/gama_gold_ddp.fits'

dat = Table.read(fpath)
# dat = dat[:30000]

assert 'DDP1_DENS' in dat.meta

# add cartesian coordinates to catalogue
xyz = cartesian(dat['RA'], dat['DEC'], dat['ZGAMA'])

dat['CARTESIAN_X'] = xyz[:,0]
dat['CARTESIAN_Y'] = xyz[:,1]
dat['CARTESIAN_Z'] = xyz[:,2]

points      = np.c_[dat['CARTESIAN_X'], dat['CARTESIAN_Y'], dat['CARTESIAN_Z']]
points      = np.array(points, copy=True)

kd_tree_all  = KDTree(points)

# randoms.
realz  = 0

rpath  = os.environ['CSCRATCH'] + '/desi/BGS/Sam/randoms_bd_{}_{:d}.fits'.format(field, realz)
rand, rand_hdr = fitsio.read(rpath, header=True)

# rand   = rand[:1000]

rpoints = np.c_[rand['CARTESIAN_X'], rand['CARTESIAN_Y'], rand['CARTESIAN_Z']]
rpoints = np.array(rpoints, copy=True)

print('Creating big rand. tree.')

big_tree = KDTree(rpoints)

print('Querying tree for closest rand.')

dd, ii   = big_tree.query([x for x in points], k=1)

dat['RANDSEP'] = dd
dat['RANDMATCH'] = rand['RANDID'][ii]
dat['BOUND_DIST'] = rand['BOUND_DIST'][ii]

nrand8 = rand_hdr['NRAND8']

dat['FILLFACTOR'] = rand['N8'][ii] / nrand8

for idx in range(3):
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

dat = dat[dat['ZGAMA'] > dat.meta['DDP1_ZMIN']]
dat = dat[dat['ZGAMA'] < dat.meta['DDP1_ZMAX']]

tiers = delta8_tier(dat['DDP1_DELTA8'])
utiers = np.unique(tiers).tolist()
utiers.remove(-99)

dat['DDP1_DELTA8_TIER'] = tiers

print(utiers)
print(np.arange(4))

assert np.all(utiers == np.arange(4))

print('Delta8 spans {} to {} over {} tiers.'.format(dat['DDP1_DELTA8'].min(), dat['DDP1_DELTA8'].max(), utiers))

for tier in utiers:
    # E.g. /global/cscratch1/sd/mjwilson/norberg//GAMA4/gama_gold_G9_ddp_n8_d0_0.fits
    isin  = (tiers == tier)
    
    opath = fpath.replace('ddp', 'ddp_n8_d0_{:d}'.format(tier))

    print('Writing {}.'.format(opath))
    
    to_write = dat[isin]
    to_write.write(opath, format='fits', overwrite=True)
    
    to_write['FIELD'] = gama_field(to_write['RA'], to_write['DEC'])
    
    isin = to_write['FIELD'] == field
    to_write_field = to_write[isin]
    
    opath_field = opath.replace('gold', 'gold_{}'.format(field))
    
    to_write_field.write(opath_field, format='fits', overwrite=True)
