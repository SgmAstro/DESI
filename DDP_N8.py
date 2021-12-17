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

# add cartesian coordinates to catalogue

xyz = cartesian(dat['RA'], dat['DEC'], dat['ZGAMA'])

dat['CARTESIAN_X'] = xyz[:,0]
dat['CARTESIAN_Y'] = xyz[:,1]
dat['CARTESIAN_Z'] = xyz[:,2]


points = np.c_[dat['CARTESIAN_X'], dat['CARTESIAN_Y'], dat['CARTESIAN_Z']]
kd_tree_all  = KDTree(points)

DDP1 = dat[dat['DDP'][:,0] == 1]
DDP2 = dat[dat['DDP'][:,1] == 1]
DDP3 = dat[dat['DDP'][:,2] == 1]

points_DDP1 = np.c_[DDP1['CARTESIAN_X'], DDP1['CARTESIAN_Y'], DDP1['CARTESIAN_Z']]
points_DDP2 = np.c_[DDP2['CARTESIAN_X'], DDP2['CARTESIAN_Y'], DDP2['CARTESIAN_Z']]
points_DDP3 = np.c_[DDP3['CARTESIAN_X'], DDP3['CARTESIAN_Y'], DDP3['CARTESIAN_Z']]

kd_tree_DDP1  = KDTree(points_DDP1)
kd_tree_DDP2  = KDTree(points_DDP2)
kd_tree_DDP3  = KDTree(points_DDP3)

indexes_DDP1  = kd_tree_all.query_ball_tree(kd_tree_DDP1, r=8.)
#indexes_DDP2  = kd_tree_all.query_ball_tree(kd_tree_DDP2, r=8.)
#indexes_DDP3  = kd_tree_all.query_ball_tree(kd_tree_DDP3, r=8.)

dat['N8_DDP1'] = [len(idx) for idx in indexes_DDP1]
#dat['N8_DDP2'] = [len(idx) for idx in indexes_DDP2]
#dat['N8_DDP3'] = [len(idx) for idx in indexes_DDP3]