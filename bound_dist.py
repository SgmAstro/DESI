import os
import gc
import sys
import tqdm
import time
import fitsio
import argparse
import multiprocessing
import numpy             as np
import astropy.io.fits   as fits
import matplotlib.pyplot as plt

from   scipy.spatial   import KDTree
from   astropy.table   import Table
from   multiprocessing import Pool
from   runtime         import calc_runtime

'''
Script to calculate the maximum distance [Mpc/h] of each random from the boundary. 
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

start  = time.time()

# https://www.dur.ac.uk/icc/cosma/cosma5/
nproc  = 16
realz  = 0

fpath  = os.environ['RANDOMS_DIR'] + '/{}_N8_{}_{:d}.fits'.format(prefix, field, realz)

if dryrun:
    fpath= fpath.replace('.fits', '_dryrun.fits')

opath = fpath.replace('{}_N8'.format(prefix), '{}_bd'.format(prefix))
    
if args.nooverwrite:
    if os.path.isfile(opath):
        print('{} found on disk and overwrite forbidden (--nooverwrite).'.format(opath))
        exit(0)

# Output is sorted by fillfactor.py;   
rand      = Table.read(fpath)

runtime   = calc_runtime(start, 'Reading {:.2f}M randoms'.format(len(rand) / 1.e6), xx=rand)
    
body      = rand[rand['IS_BOUNDARY'] == 0]
boundary  = rand[rand['IS_BOUNDARY'] == 1]

split_idx = np.arange(len(body))
splits    = np.array_split(split_idx, 10 * nproc)

bids      = boundary['RANDID']

boundary  = np.c_[boundary['CARTESIAN_X'], boundary['CARTESIAN_Y'], boundary['CARTESIAN_Z']]
boundary  = np.array(boundary, copy=True)

kd_tree   = KDTree(boundary)

runtime   = calc_runtime(start, 'Created boundary tree.')

# points  = np.c_[body['CARTESIAN_X'], body['CARTESIAN_Y'], body['CARTESIAN_Z']] 
# points  = [x for x in points]
# dd, ii  = kd_tree.query(points, k=1)

del rand
del boundary
del split_idx

gc.collect()

'''
local_vars = list(locals().items())                                                                                                                                                                  

for var, obj in local_vars:                                                                                                                                                                           
    print(var, sys.getsizeof(obj))                                                                                                                                                                   

exit()
'''

runtime   = calc_runtime(start, 'Deleted rand.')

def process_one(split, pid=0):
    points = np.c_[body[split]['CARTESIAN_X'], body[split]['CARTESIAN_Y'], body[split]['CARTESIAN_Z']] 
    points = [x for x in points]

    try:
        pid  = multiprocessing.current_process().name.ljust(20)

    except Exception as e:
        print(e)


    # print('Querying boundary tree for split [{}...{}]'.format(split[0], split[-1]))
    
    dd, ii = kd_tree.query(points, k=1)
    
    del  points

    return  dd.tolist(), ii.tolist()

'''
# Serial.

result = []

for split in splits:
    result += process_one(split)

result = np.array(result)
'''

runtime = calc_runtime(start, 'POOL:  Querying bound dist for body points.')

with Pool(nproc) as pool:
    # result  = p.map(process_one, splits)

    results = []

    for result in tqdm.tqdm(pool.imap(process_one, iterable=splits), total=len(splits)):
        results.append(result)

    pool.close()
    pool.join()

runtime = calc_runtime(start, 'POOL:  Done with queries')

flat_result = []
flat_ii     = []

for rr in results:
    flat_result += rr[0]
    flat_ii     += rr[1]

rand = Table.read(fpath)
rand['BOUND_DIST'] = 0.0
rand['BOUND_ID']   = 0

rand['BOUND_DIST'][rand['IS_BOUNDARY'] == 0] = np.array(flat_result)
rand['BOUND_ID'][rand['IS_BOUNDARY'] == 0]   = bids[np.array(flat_ii)]

sphere_radius = rand.meta['RSPHERE']
rand['FILLFACTOR']   = np.clip(rand['FILLFACTOR'], 0., 1.)
rand['FILLFACTOR'][rand['BOUND_DIST'].data > sphere_radius] = 1

runtime = calc_runtime(start, 'Shuffling')

# randomise rows.                                                                                                                                                
idx  = np.arange(len(rand))
idx  =  np.random.choice(idx, size=len(idx), replace=False)

rand = rand[idx]

# Bound dist.
# https://docs.scipy.org/doc/scipy/reference/generated/scipy.spatial.KDTree.query.html#scipy.spatial.KDTree.query

runtime = calc_runtime(start, 'Writing {}'.format(opath), xx=rand)

rand.write(opath, format='fits', overwrite=True)

runtime   = calc_runtime(start, 'Finished')


