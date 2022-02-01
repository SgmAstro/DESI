import os
import gc
import sys
import time
import tqdm
import fitsio
import argparse
import multiprocessing
import numpy               as np
import astropy.io.fits     as fits
import matplotlib.pyplot   as plt

from   scipy.spatial       import KDTree
from   astropy.table       import Table
from   multiprocessing     import Pool
from   runtime             import calc_runtime
from   memory_profiler     import profile


parser = argparse.ArgumentParser(description='Calculate fill factor using randoms.')
parser.add_argument('-f', '--field', type=str, help='Sselect equatorial GAMA field: G9, G12, G15', default='G9')
parser.add_argument('-d', '--dryrun', help='Dryrun.', action='store_true')
parser.add_argument('--prefix', help='filename prefix', default='randoms')
parser.add_argument('--nproc', help='nproc', default=16, type=np.int32)
parser.add_argument('--maxtasksperchild', help='maxtasksperchild', default=1000, type=np.int32)
parser.add_argument('--realz', help='Realization number', default=0, type=np.int32)
parser.add_argument('--nooverwrite',  help='Do not overwrite outputs if on disk', action='store_true')

args   = parser.parse_args()

field  = args.field.upper()
dryrun = args.dryrun
prefix = args.prefix

maxtasksperchild = args.maxtasksperchild

# https://www.dur.ac.uk/icc/cosma/cosma5/
nproc  = args.nproc
realz  = args.realz

fpath  = os.environ['RANDOMS_DIR'] + '/{}_{}_{:d}.fits'.format(prefix, field, realz)
start  = time.time()

if dryrun:
    fpath = fpath.replace('.fits', '_dryrun.fits')

opath  = fpath.replace('{}_{}'.format(prefix, field), '{}_N8_{}'.format(prefix, field))

if args.nooverwrite:
    if os.path.isfile(fpath) and os.path.isfile(opath):
        print('{} found on disk and overwrite forbidden (--nooverwrite).'.format(fpath))
        print('{} found on disk and overwrite forbidden (--nooverwrite).'.format(opath))
        exit(0)
    
# Read randoms file, split by field (DDP1, or not).
rand      = Table.read(fpath)
    
runtime   = calc_runtime(start, 'Reading {:.2f}M randoms'.format(len(rand) / 1.e6), xx=rand)

rand.sort('CARTESIAN_X')

runtime   = calc_runtime(start, 'Sorted randoms by X')

split_idx = np.arange(len(rand))
splits    = np.array_split(split_idx, 10 * nproc)

runtime   = calc_runtime(start, 'Split randoms by {} nproc'.format(nproc))

points    = np.c_[rand['CARTESIAN_X'], rand['CARTESIAN_Y'], rand['CARTESIAN_Z']]

runtime   = calc_runtime(start, 'Creating big tree.')

big_tree  = KDTree(points)
runtime   = calc_runtime(start, 'Created big (randoms) tree')

del rand
del points
del split_idx

gc.collect()

runtime   = calc_runtime(start, 'Deleted rand.')

'''
local_vars = list(locals().items())

for var, obj in local_vars:
    print(var, sys.getsizeof(obj))
'''
    
def process_one(split, pid=0):
    _points  = np.c_[big_tree.data[split,0], big_tree.data[split,1], big_tree.data[split,2]] 
    _points  = np.array(_points, copy=True)

    try:
        pid  = multiprocessing.current_process().name.ljust(20)

    except Exception as e:
        print(e)
    
    # msg      = 'POOL {}:  Creating split [{} ... {}] tree.'.format(pid, split[0], split[-1])
    # runtime  = calc_runtime(start, msg)
        
    kd_tree  = KDTree(_points)

    # msg      = 'POOL {}:  Querying split [{} ... {}] tree.'.format(pid, split[0], split[-1])
    # runtime  = calc_runtime(start, msg)

    # https://docs.scipy.org/doc/scipy/reference/generated/scipy.spatial.KDTree.query_ball_tree.html#scipy.spatial.KDTree.query_ball_tree
    # https://docs.scipy.org/doc/scipy/reference/generated/scipy.spatial.KDTree.count_neighbors.html#scipy.spatial.KDTree.count_neighbors
    indexes  = kd_tree.query_ball_tree(big_tree, r=8.)

    del  _points
    del  kd_tree
   
    # runtime  = calc_runtime(start, 'Flattening')

    flat     = [len(idx) for idx in indexes]

    # runtime  = calc_runtime(start, 'Flattened')

    return  flat


runtime = calc_runtime(start, 'POOL:  Counting < 8 Mpc/h pairs for small trees.')

# maxtasksperchild=maxtasksperchild
with Pool(nproc) as pool:
    # result = pool.map(process_one,  splits)
    # result = pool.imap(process_one, splits)

    results = []

    for result in tqdm.tqdm(pool.imap(process_one, iterable=splits), total=len(splits)):
        results.append(result)

    pool.close()
    pool.join()

    '''
    local_vars = list(locals().items())                                                                                                                                                            
    
    for var, obj in local_vars:                                                                                                                                                                     
        print(var, sys.getsizeof(obj)) 
    '''

runtime = calc_runtime(start, 'POOL:  Done with queries')

flat_result = []
    
for rr in results:
    flat_result += rr

rand                 = Table.read(fpath)
rand.sort('CARTESIAN_X')

rand['RAND_N8']      = np.array(flat_result).astype(np.int32)
rand['FILLFACTOR']   = rand['RAND_N8'] / rand.meta['NRAND8']

rand.meta['RSPHERE'] = 8.
    
# TODO: INHERIT FILL FACTOR THRESHOLD FROM PARAMS FILE.
rand.meta['FILLFACTOR_INFRAC'] = np.mean(rand['FILLFACTOR'] > 0.8)

runtime = calc_runtime(start, 'Writing {}.'.format(opath), xx=rand)

rand.write(opath, format='fits', overwrite=True)

runtime = calc_runtime(start, 'Finished')
