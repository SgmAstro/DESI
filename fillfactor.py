import os
import fitsio
import numpy as np
import astropy.io.fits as fits
import matplotlib.pyplot as plt

from   scipy.spatial import KDTree
from   astropy.table import Table
from   multiprocessing import Pool

import argparse


parser = argparse.ArgumentParser(description='Select GAMA field.')
parser.add_argument('-f', '--field', type=str, help='select equatorial GAMA field: G9, G12, G15', required=True)
parser.add_argument('-d', '--dryrun', help='Dryrun.', action='store_true')

args = parser.parse_args()

field = args.field.upper()
dryrun = args.dryrun

nproc = 16
realz = 0

fpath = os.environ['RANDOMS_DIR'] + '/randoms_{}_{:d}.fits'.format(field, realz)

print('Reading rand.')

if dryrun:
    fpath = fpath.replace('.fits', '_dryrun.fits')

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

rand['N8'] = np.array(flat_result).astype(np.int32)
rand.meta['RSPHERE'] = 8.

opath = fpath.replace('randoms_{}'.format(field), 'randoms_N8_{}'.format(field))

print('Writing {}.'.format(opath))

rand.write(opath, format='fits', overwrite=True)
