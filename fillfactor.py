import os
import fitsio
import numpy as np
import astropy.io.fits as fits
import matplotlib.pyplot as plt

from   scipy.spatial import KDTree
from   astropy.table import Table
from   multiprocessing import Pool

nproc = 8
fpath = os.environ['CSCRATCH'] + '/desi/BGS/Sam/randoms.fits'

rand  = fits.open(fpath)
hdr   = rand[1].header

# print(rand.info())
# print(hdr)

rand     = fitsio.read(fpath)
# rand = rand[:200*nproc]

split_idx = np.arange(len(rand))
splits = np.split(split_idx, nproc)

def process_one(split):
    sub_rand = rand[split]
    sub_rand = Table(sub_rand, copy=True)

    points   = np.c_[sub_rand['CARTESIAN_X'], sub_rand['CARTESIAN_Y'], sub_rand['CARTESIAN_Z']] 

    kd_tree  = KDTree(points)

    indexes  = kd_tree.query_ball_tree(kd_tree, r=8.)

    del sub_rand
    del kd_tree
    
    return  [len(idx) for idx in indexes]

#print(splits)
#print(rand[splits[0]])
'''
result = []

for split in splits:
    result.append(process_one(split))

result = np.array(result)
'''
#print(result)

with Pool(nproc) as p:
        result = p.map(process_one, splits)

#print(result)

flat_result = []

for rr in result:
    flat_result += rr

rand = Table(rand)
rand['N8'] = np.array(flat_result)

# Bound dist.
# https://docs.scipy.org/doc/scipy/reference/generated/scipy.spatial.KDTree.query.html#scipy.spatial.KDTree.query

rand.write(fpath.replace('randoms', 'randoms_N8'), format='fits', overwrite=True)
