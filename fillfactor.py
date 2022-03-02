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

from   findfile            import findfile, fetch_fields, overwrite_check

parser = argparse.ArgumentParser(description='Calculate fill factor using randoms.')
parser.add_argument('-f', '--field', type=str, help='Select equatorial GAMA field: G9, G12, G15', default='G9')
parser.add_argument('-d', '--dryrun', help='Dryrun.', action='store_true')
parser.add_argument('-s', '--survey', help='Select survey.', default='gama')
parser.add_argument('--prefix', help='filename prefix', default='randoms')
parser.add_argument('--nproc', help='nproc', default=8, type=int)
parser.add_argument('--maxtasksperchild', help='maxtasksperchild', default=1000, type=np.int32)
parser.add_argument('--realz', help='Realization number', default=0, type=np.int32)
parser.add_argument('--nooverwrite',  help='Do not overwrite outputs if on disk', action='store_true')
parser.add_argument('--max_nsplit',  default=-99, type=int)

args      x = parser.parse_args()

field      = args.field.upper()
dryrun     = args.dryrun
prefix     = args.prefix
survey     = args.survey.lower()
max_nsplit = args.max_nsplit

maxtasksperchild = args.maxtasksperchild

fields = fetch_fields(survey)

assert field in fields, 'Error: Field not in fields'

# https://www.dur.ac.uk/icc/cosma/cosma5/
nproc  = args.nproc
realz  = args.realz

fpath  = findfile(ftype='randoms', dryrun=dryrun, field=field, survey=survey, prefix=prefix)

start  = time.time()

opath  = findfile(ftype='randoms_n8', dryrun=dryrun, field=field, survey=survey, prefix=prefix)

if args.nooverwrite:
    overwrite_check(opath)
    
# Read randoms file, split by field (DDP1, or not).
rand      = Table.read(fpath)
    
runtime   = calc_runtime(start, 'Reading {:.2f}M randoms'.format(len(rand) / 1.e6), xx=rand)

rand.sort('CARTESIAN_X')

runtime   = calc_runtime(start, 'Sorted randoms by X')

splits    = np.arange(len(rand))

if max_nsplit > 0:
    splits = splits[:max_nsplit] 
    
    print('Warning:  limiting splits to {}'.format(len(splits)))

runtime   = calc_runtime(start, 'Split randoms by {} nproc'.format(nproc))

points    = np.c_[rand['CARTESIAN_X'].data.astype(np.float32), rand['CARTESIAN_Y'].data.astype(np.float32), rand['CARTESIAN_Z'].data.astype(np.float32)]

runtime   = calc_runtime(start, 'Creating big tree.')

big_tree  = KDTree(points, leafsize=5)
runtime   = calc_runtime(start, 'Created big (randoms) tree')

del rand
del points
del split_idx

gc.collect()

runtime   = calc_runtime(start, 'Deleted rand.')

def process_one(split, pid=0):
    _points    = np.c_[big_tree.data[split,0], big_tree.data[split,1], big_tree.data[split,2]] 

    # HACK
    # _points  = np.array(_points, copy=True)

    try:
        pid    = multiprocessing.current_process().name.ljust(20)

    except Exception as e:
        print(e)
    
    msg      = 'POOL {}:  Creating {} long split [{} ... {}] tree.'.format(pid, len(split), split[0], split[-1])
    runtime  = calc_runtime(start, msg)
        
    kd_tree  = KDTree(_points, leafsize=5)

    msg      = 'POOL {}:  Querying {} long split [{} ... {}] tree.'.format(pid, len(split), split[0], split[-1])
    runtime  = calc_runtime(start, msg)

    # https://docs.scipy.org/doc/scipy/reference/generated/scipy.spatial.KDTree.query_ball_tree.html#scipy.spatial.KDTree.query_ball_tree
    # https://docs.scipy.org/doc/scipy/reference/generated/scipy.spatial.KDTree.count_neighbors.html#scipy.spatial.KDTree.count_neighbors
    indexes  = kd_tree.query_ball_tree(big_tree, r=8.)

    del  _points
    del  kd_tree
   
    # runtime  = calc_runtime(start, 'Flattening')

    flat     = [len(idx) for idx in indexes]

    # runtime  = calc_runtime(start, 'Flattened')

    return  flat

chunksize   = 1000
nchunk      = 1. * len(splits) / chunksize

runtime     = calc_runtime(start, 'POOL:  Counting < 8 Mpc/h pairs for small trees of {} splits.'.format(nchunk))

pool_start  = time.time()

results     = [process_one(splits[0], pid=0)]

split_time  = time.time() - pool_start
split_time /= 60.

runtime     = calc_runtime(start, 'POOL:  Expected runtime of {:.6f} minutes with {:d} proc. and split time {:.6f} mins'.format(nchunk * split_time / nproc, nproc, split_time))

done_nsplit = 1

# maxtasksperchild=maxtasksperchild
with Pool(nproc) as pool:
    # result = pool.map(process_one,  splits)
    # result = pool.imap(process_one, splits)

    for result in tqdm.tqdm(pool.imap(process_one, iterable=splits[1:], chunksize=chunksize), total=(nchunk-1))):
        results.append(result)

        done_nsplit  += 1
        
        if (done_nsplit % nproc) == 0:
            pool_time = (time.time() - pool_start)  / 60.
            runtime   = calc_runtime(start, 'POOL:  New expected runtime of {:.3f} minutes with {:d} proc.'.format(len(splits) * pool_time / done_nsplit, nproc))

    pool.close()
    pool.join()

runtime     = calc_runtime(start, 'POOL:  Done with queries of {} splits with effective split time {}'.format(done_nsplit, pool_time / done_nsplit))

flat_result = []
    
for rr in results:
    flat_result += rr

rand                 = Table.read(fpath)
rand.sort('CARTESIAN_X')

rand['RAND_N8']      = np.array(flat_result).astype(np.int32)
rand['FILLFACTOR']   = rand['RAND_N8'] / rand.meta['NRAND8']
rand.meta['RSPHERE'] = 8.

boundary = Table.read(fpath, 'BOUNDARY')

header   = fits.Header()

hx       = fits.HDUList()
hx.append(fits.PrimaryHDU(header=header))
hx.append(fits.convenience.table_to_hdu(rand))
hx.append(fits.convenience.table_to_hdu(boundary))

runtime = calc_runtime(start, 'Writing {}.'.format(opath), xx=rand)

hx.writeto(opath, overwrite=True)

runtime = calc_runtime(start, 'Finished')
