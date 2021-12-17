import os
import fitsio
import numpy as np
import astropy.io.fits as fits
import matplotlib.pyplot as plt

from   scipy.spatial import KDTree
from   astropy.table import Table
from   multiprocessing import Pool


nproc = 16

field = 'G9'
realz = 0

dryrun=False

fpath = os.environ['CSCRATCH'] + '/desi/BGS/Sam/randoms_{}_{:d}.fits'.format(field, realz)

print('Reading rand.')

rand = Table.read(fpath)

if dryrun:
    rand = rand[:200*nproc]
    fpath = fpath.replace('.fits', '_dryrun.fits')
    
sort_idx = np.argsort(rand['CARTESIAN_X'])
rand = rand[sort_idx]

split_idx = np.arange(len(rand))
splits = np.array_split(split_idx, nproc)

points = np.c_[rand['CARTESIAN_X'], rand['CARTESIAN_Y'], rand['CARTESIAN_Z']]
points = np.array(points, copy=True)

print('Creating big tree.')

big_tree = KDTree(points)

def process_one(split):
    _points   = np.c_[rand[split]['CARTESIAN_X'], rand[split]['CARTESIAN_Y'], rand[split]['CARTESIAN_Z']] 
    _points = np.array(_points, copy=True)
    
    print('Creating split [{} ... {}] tree.'.format(split[0], split[-1]))
    
    kd_tree  = KDTree(_points)

    print('Querying split [{} ... {}] tree.'.format(split[0], split[-1]))
    
    indexes  = kd_tree.query_ball_tree(big_tree, r=8.)

    del kd_tree
    
    return  [len(idx) for idx in indexes]

'''
## serial.
result = []

for split in splits:
    result.append(process_one(split))

result = np.array(result)
'''
print('Counting < 8 Mpc/h pairs for small trees.')

with Pool(nproc) as p:
    result = p.map(process_one, splits)

print('Flattening.')

flat_result = []

for rr in result:
    flat_result += rr

rand['N8'] = np.array(flat_result).astype(np.int32)

print('Writing.')

rand.write(fpath.replace('randoms', 'randoms_N8'), format='fits', overwrite=True)
