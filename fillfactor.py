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


parser = argparse.ArgumentParser(description='Calculate fill factor using randoms.')
parser.add_argument('-f', '--field', type=str, help='Sselect equatorial GAMA field: G9, G12, G15', required=True)
parser.add_argument('-d', '--dryrun', help='Dryrun.', action='store_true')
parser.add_argument('--prefix', help='filename prefix', default='randoms')
parser.add_argument('--nooverwrite',  help='Do not overwrite outputs if on disk', action='store_true')

args   = parser.parse_args()

field  = args.field.upper()
dryrun = args.dryrun
prefix = args.prefix

nproc = 16
realz = 0

fpath = os.environ['RANDOMS_DIR'] + '/{}_{}_{:d}.fits'.format(prefix, field, realz)

start = time.time()

print('Reading rand.')

if dryrun:
    fpath = fpath.replace('.fits', '_dryrun.fits')

opath = fpath.replace('{}_{}'.format(prefix, field), '{}_N8_{}'.format(prefix, field))

if args.nooverwrite:
    if os.path.isfile(fpath) and os.path.isfile(opath):
        print('{} found on disk and overwrite forbidden (--nooverwrite).'.format(fpath))
        print('{} found on disk and overwrite forbidden (--nooverwrite).'.format(opath))
        exit(0)
    
rand = Table.read(fpath)
    
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

print('Counting < 8 Mpc/h pairs for small trees.')

with Pool(nproc) as p:
    result = p.map(process_one, splits)

print('Flattening.')

flat_result = []

for rr in result:
    flat_result += rr

rand['RAND_N8'] = np.array(flat_result).astype(np.int32)
rand['FILLFACTOR']  = rand['RAND_N8'] / rand.meta['NRAND8']

# TODO:  DDP1_FILLFACTOR.

rand.meta['RSPHERE'] = 8.

# TODO: INHERIT FILL FACTOR THRESHOLD FROM PARAMS FILE.
rand.meta['FILLFACTOR_INFRAC'] = np.mean(rand['FILLFACTOR'] > 0.8)

print('Writing {}.'.format(opath))

rand.write(opath, format='fits', overwrite=True)

print('Finished in {} mins.'.format((time.time() - start) / 60.))
