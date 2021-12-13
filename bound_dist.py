import os
import fitsio
import numpy as np
import astropy.io.fits as fits
import matplotlib.pyplot as plt

from   scipy.spatial import KDTree
from   astropy.table import Table
from   multiprocessing import Pool

nproc = 8
fpath = os.environ['CSCRATCH'] + '/desi/BGS/Sam/randoms_N8.fits'

rand  = fits.open(fpath)
hdr   = rand[1].header

# print(rand.info())
# print(hdr)

rand     = fitsio.read(fpath)
# rand = rand[:200*nproc]

body = rand[rand['IS_BOUNDARY'] == 0]
boundary = rand[rand['IS_BOUNDARY'] == 1]

split_idx = np.arange(len(body))
splits = np.array_split(split_idx, nproc)

boundary = np.c_[boundary['CARTESIAN_X'], boundary['CARTESIAN_Y'], boundary['CARTESIAN_Z']]

kd_tree  = KDTree(boundary)

def process_one(split):
    sub_rand = body[split]
    sub_rand = Table(sub_rand, copy=True)

    points   = np.c_[sub_rand['CARTESIAN_X'], sub_rand['CARTESIAN_Y'], sub_rand['CARTESIAN_Z']] 

    points = [x for x in points]

    dd, ii = kd_tree.query(points, k=1)

    # ii = boundary['RANDID'][ii]
    
    return  dd.tolist() #, ii.tolist()
    
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
    flat_result += rr
    # flat_ii += rr[1]
    
rand = Table(rand)
rand['BOUND_DIST'] = 0.0
rand['BOUND_ID'] = 0

rand['BOUND_DIST'][rand['IS_BOUNDARY'] == 0] = np.array(flat_result)
# rand['BOUND_ID'][rand['IS_BOUNDARY'] == 0] = np.array(flat_ii)

# Bound dist.
# https://docs.scipy.org/doc/scipy/reference/generated/scipy.spatial.KDTree.query.html#scipy.spatial.KDTree.query

rand.write(fpath.replace('randoms_N8', 'randoms_bd'), format='fits', overwrite=True)
