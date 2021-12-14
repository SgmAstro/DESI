import os
import fitsio
import numpy as np
import astropy.io.fits as fits
import matplotlib.pyplot as plt

from   scipy.spatial import KDTree
from   astropy.table import Table
from   multiprocessing import Pool

np.random.seed(314)

nproc = 12

field = 'G9'
realz = 0

fpath = os.environ['CSCRATCH'] + '/desi/BGS/Sam/randoms_N8_{}_{:d}.fits'.format(field, realz)

#rand  = fits.open(fpath)
#rand  =rand[:200]

# hdr   = rand[1].header

# print(rand.info())
# print(hdr)

# Outpute is sorted by fillfactor.py;
rand     = fitsio.read(fpath)
#rand = rand[:800*nproc]

body = rand[rand['IS_BOUNDARY'] == 0]
boundary = rand[rand['IS_BOUNDARY'] == 1]

split_idx = np.arange(len(body))
splits = np.array_split(split_idx, nproc)

bids = boundary['RANDID']

boundary = np.c_[boundary['CARTESIAN_X'], boundary['CARTESIAN_Y'], boundary['CARTESIAN_Z']]

print('Creating boundary tree')

kd_tree  = KDTree(boundary)

points = np.c_[body['CARTESIAN_X'], body['CARTESIAN_Y'], body['CARTESIAN_Z']] 
points = [x for x in points]

# dd, ii = kd_tree.query(points, k=1)

def process_one(split):
    # sub_rand = Table(body[split], copy=True)

    points = np.c_[body[split]['CARTESIAN_X'], body[split]['CARTESIAN_Y'], body[split]['CARTESIAN_Z']] 
    points = [x for x in points]

    print('Querying boundary tree for split [{}...{}]'.format(split[0], split[-1]))
    
    dd, ii = kd_tree.query(points, k=1)
    
    return  dd.tolist(), ii.tolist()

#print(splits)
#print(rand[splits[0]])
'''
result = []

for split in splits:
    result += process_one(split)

result = np.array(result)
'''
#print(result)

with Pool(nproc) as p:
    result = p.map(process_one, splits)

#print(result)

flat_result = []
flat_ii = []

for rr in result:
    flat_result += rr[0]
    flat_ii += rr[1]

rand = Table(rand)
rand['BOUND_DIST'] = 0.0
rand['BOUND_ID'] = 0

rand['BOUND_DIST'][rand['IS_BOUNDARY'] == 0] = np.array(flat_result)
rand['BOUND_ID'][rand['IS_BOUNDARY'] == 0] = bids[np.array(flat_ii)]

print('Shuffling.')

# randomise rows.                                                                                                                                                
idx = np.arange(len(rand))
idx=  np.random.choice(idx, size=len(idx), replace=False)

rand = rand[idx]

# Bound dist.
# https://docs.scipy.org/doc/scipy/reference/generated/scipy.spatial.KDTree.query.html#scipy.spatial.KDTree.query

print('Writing {}.'.format(fpath.replace('randoms_N8', 'randoms_bd')))

rand.write(fpath.replace('randoms_N8', 'randoms_bd'), format='fits', overwrite=True)
