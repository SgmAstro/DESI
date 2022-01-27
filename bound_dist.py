import os
import time
import fitsio
import numpy as np
import astropy.io.fits as fits
import matplotlib.pyplot as plt

from   scipy.spatial import KDTree
from   astropy.table import Table
from   multiprocessing import Pool

import argparse

'''
Script to calculate the maximimum distance [Mpc/h] of each random from the boundary. 
'''

np.random.seed(314)

parser = argparse.ArgumentParser(description='Find boundary distance for all randoms in a specified field..')
parser.add_argument('-f', '--field', type=str, help='Select equatorial GAMA field: G9, G12, G15', required=True)
parser.add_argument('-d', '--dryrun', help='Dryrun.', action='store_true')
parser.add_argument('--prefix', help='filename prefix', default='randoms')
parser.add_argument('--nooverwrite',  help='Do not overwrite outputs if on disk', action='store_true')

args   = parser.parse_args()

field  = args.field.upper()
dryrun = args.dryrun
prefix = args.prefix

start = time.time()

nproc = 12
realz = 0

fpath = os.environ['RANDOMS_DIR'] + '/{}_N8_{}_{:d}.fits'.format(prefix, field, realz)

if dryrun:
    fpath= fpath.replace('.fits', '_dryrun.fits')

opath = fpath.replace('{}_N8'.format(prefix), '{}_bd'.format(prefix))
    
if args.nooverwrite:
    if os.path.isfile(fpath) and os.path.isfile(opath):
        print('{} found on disk and overwrite forbidden (--nooverwrite).'.format(fpath))
        print('{} found on disk and overwrite forbidden (--nooverwrite).'.format(opath))

        exit(0)

# Output is sorted by fillfactor.py;   
rand  = Table.read(fpath)
    
body = rand[rand['IS_BOUNDARY'] == 0]
boundary = rand[rand['IS_BOUNDARY'] == 1]

split_idx = np.arange(len(body))
splits = np.array_split(split_idx, nproc)

bids = boundary['RANDID']

boundary = np.c_[boundary['CARTESIAN_X'], boundary['CARTESIAN_Y'], boundary['CARTESIAN_Z']]
boundary = np.array(boundary, copy=True)

print('Creating boundary tree')

kd_tree  = KDTree(boundary)

# points = np.c_[body['CARTESIAN_X'], body['CARTESIAN_Y'], body['CARTESIAN_Z']] 
# points = [x for x in points]

# dd, ii = kd_tree.query(points, k=1)

def process_one(split):
    points = np.c_[body[split]['CARTESIAN_X'], body[split]['CARTESIAN_Y'], body[split]['CARTESIAN_Z']] 
    points = [x for x in points]

    print('Querying boundary tree for split [{}...{}]'.format(split[0], split[-1]))
    
    dd, ii = kd_tree.query(points, k=1)
    
    return  dd.tolist(), ii.tolist()

'''
# Serial.

result = []

for split in splits:
    result += process_one(split)

result = np.array(result)
'''
with Pool(nproc) as p:
    result = p.map(process_one, splits)

flat_result = []
flat_ii = []

for rr in result:
    flat_result += rr[0]
    flat_ii += rr[1]

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

print('Writing {}.'.format(opath))

rand.write(opath, format='fits', overwrite=True)

print('Finished in {} mins.'.format((time.time() - start) / 60.))
