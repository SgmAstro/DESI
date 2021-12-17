import os
import fitsio
import numpy as np
import matplotlib.pyplot as plt

from astropy.table import Table
from scipy.spatial import KDTree
from cartesian import cartesian

fpath = os.environ['CSCRATCH'] + '/norberg/GAMA4/gama_gold_ddp.fits'

dat = fitsio.read(fpath)
dat = Table(dat)

dat = dat[:1000]

# add cartesian coordinates to catalogue
xyz = cartesian(dat['RA'], dat['DEC'], dat['ZGAMA'])

dat['CARTESIAN_X'] = xyz[:,0]
dat['CARTESIAN_Y'] = xyz[:,1]
dat['CARTESIAN_Z'] = xyz[:,2]

points       = np.c_[dat['CARTESIAN_X'], dat['CARTESIAN_Y'], dat['CARTESIAN_Z']]
points       = np.array(points, copy=True)

kd_tree_all  = KDTree(points)

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
