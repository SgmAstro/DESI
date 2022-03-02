import os
import numpy as np
import matplotlib.pyplot as plt
import fitsio
from astropy.table import Table
import cosmo as cosmo
from   cartesian import cartesian

from   scipy.spatial import KDTree
from   astropy.table import Table
from   multiprocessing import Pool

nproc = 10
field = 'G9'
realz = 0
boundary_percent = 1.
area = 12*5

fpath = os.environ['CSCRATCH'] + '/desi/BGS/Sam/randoms_{}_{:d}.fits'.format(field, realz)
rand = fitsio.read(fpath)
rand = rand[:200*nproc]
sort_idx = np.argsort(rand['CARTESIAN_X'])


# load in the GAMA galaxies
fpath = os.environ['CSCRATCH'] + '/norberg/GAMA4/gama_gold_kE_5k.fits'.format(field, realz)
dat = fitsio.read(fpath)

# set to G9 field for now
dat = dat[dat['RA'] < 160]
ras  = dat['RA']
decs = dat['DEC']
zs   = dat['ZGAMA']
xyz  = cartesian(ras, decs, zs)

dat = Table(dat)

dat['CARTESIAN_X'] = xyz[:,0]
dat['CARTESIAN_Y'] = xyz[:,1]
dat['CARTESIAN_Z'] = xyz[:,2]
dat['V'] = cosmo.volcom(dat['ZGAMA'], area)

dat['IS_BOUNDARY'] = 0

dat['IS_BOUNDARY'][dat['RA']  > np.percentile(dat['RA'], 100. - boundary_percent)] = 1
dat['IS_BOUNDARY'][dat['RA']  < np.percentile(dat['RA'],  boundary_percent)] = 1

dat['IS_BOUNDARY'][dat['DEC'] > np.percentile(dat['DEC'], 100. - boundary_percent)] = 1
dat['IS_BOUNDARY'][dat['DEC'] < np.percentile(dat['DEC'], boundary_percent)] = 1

dat['IS_BOUNDARY'][dat['V'] >= np.percentile(dat['V'], 100. - boundary_percent)] = 1
dat['IS_BOUNDARY'][dat['V'] <= np.percentile(dat['V'],  boundary_percent)] = 1

points_gama = np.c_[dat['CARTESIAN_X'], dat['CARTESIAN_Y'], dat['CARTESIAN_Z']]
points_rand = np.c_[rand['CARTESIAN_X'], rand['CARTESIAN_Y'], rand['CARTESIAN_Z']]

# calculate bound_dist here
body = dat[dat['IS_BOUNDARY'] == 0]
boundary = dat[dat['IS_BOUNDARY'] == 1]
boundary = np.c_[boundary['CARTESIAN_X'], boundary['CARTESIAN_Y'], boundary['CARTESIAN_Z']]

kd_tree  = KDTree(boundary)

points = np.c_[body['CARTESIAN_X'], body['CARTESIAN_Y'], body['CARTESIAN_Z']] 
points = [x for x in points]

dd, ii = kd_tree.query(points, k=1)

dat = Table(dat)
dat['BOUND_DIST'] = 0.0

dat['BOUND_DIST'][dat['IS_BOUNDARY'] == 0] = np.array(dd)


plt.figure(figsize=(6, 6))
plt.plot(points_rand[:, 0], points_rand[:, 1], "xk", markersize=14)
plt.plot(points_gama[:, 0], points_gama[:, 1], "og", markersize=14)
kd_tree_rand = KDTree(points_rand)
kd_tree_gama = KDTree(points_gama)

indexes = kd_tree_gama.query_ball_tree(kd_tree_rand, r=8)
N8_dr = [len(idx) for idx in indexes]
dat['N8'] = N8_dr