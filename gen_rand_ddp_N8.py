import os
import fitsio
import numpy as np
import matplotlib.pyplot as plt

from astropy.table import Table
from scipy.spatial import KDTree
from cartesian import cartesian
from delta8_limits import dd8_limits, delta8_tier


field = 'G9'
realz = 0

fpath = os.environ['CSCRATCH'] + '/norberg/GAMA4/gama_gold_ddp_n8.fits'
dat = Table.read(fpath)
# dat = dat[:1000]

fpath = os.environ['CSCRATCH'] + '/desi/BGS/Sam/randoms_bd_{}_{:d}.fits'.format(field, realz)
rand = Table.read(fpath)
# rand = rand[:1000]

points       = np.c_[rand['CARTESIAN_X'], rand['CARTESIAN_Y'], rand['CARTESIAN_Z']]
points       = np.array(points, copy=True)

kd_tree_rand = KDTree(points)

for idx in range(3):
    ddp_idx      = idx + 1

    print('Solving for DDP {}'.format(ddp_idx))
    
    ddp          = dat[dat['DDP'][:,idx] == 1]
    points_ddp   = np.c_[ddp['CARTESIAN_X'], ddp['CARTESIAN_Y'], ddp['CARTESIAN_Z']]
    points_ddp   = np.array(points_ddp, copy=True)
    
    kd_tree_ddp  = KDTree(points_ddp)
    
    indexes_ddp  = kd_tree_rand.query_ball_tree(kd_tree_ddp, r=8.)

    rand['DDP{:d}_N8'.format(ddp_idx)] = np.array([len(idx) for idx in indexes_ddp])

# HACK: can be removed after full randoms run through.
rand.meta['NRAND8'] = 1072.33029242

rand['FILLFACTOR']  = rand['N8'] / rand.meta['NRAND8']   

rand['DDP1_DELTA8'] = (rand['DDP1_N8'] / (dat.meta['VOL8'] * dat.meta['DDP1_DENS']) / rand['FILLFACTOR']) - 1.
rand['DDP2_DELTA8'] = (rand['DDP2_N8'] / (dat.meta['VOL8'] * dat.meta['DDP2_DENS']) / rand['FILLFACTOR']) - 1.
rand['DDP3_DELTA8'] = (rand['DDP3_N8'] / (dat.meta['VOL8'] * dat.meta['DDP3_DENS']) / rand['FILLFACTOR']) - 1.

rand['DDP1_DELTA8_TIER'] = delta8_tier(rand['DDP1_DELTA8'])

print('Writing {}'.format(fpath.replace('bd', 'bd_ddp_n8')))
    
rand.write(fpath.replace('bd', 'bd_ddp_n8'), format='fits', overwrite=True)
