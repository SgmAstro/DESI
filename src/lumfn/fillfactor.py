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

from   functools           import partial
from   datetime            import datetime
from   scipy.spatial       import KDTree
from   astropy.table       import Table
from   multiprocessing     import Pool
from   runtime             import calc_runtime
from   findfile            import findfile, fetch_fields, overwrite_check, call_signature, gather_cat
from   config              import Configuration
from   ddp_zlimits         import ddp_zlimits
from   params              import sphere_radius


def collate_fillfactors(realzs=np.array([0]), field='G9', survey='gama', dryrun=False, prefix=None, write=True, force=False, oversample=2):
    print('Collating fillfactor realizations into main (realz=0).')
    
    realzs     = np.sort(realzs)

    assert realzs[0] == 0

    opaths     = [findfile(ftype='randoms_n8', dryrun=dryrun, field=field, survey=survey, prefix=prefix, realz=realz) for realz in realzs]
    opath      = opaths[0] 

    assert  np.all(np.array([os.path.isfile(x) for x in opaths])), 'Failed to find {}'.format(opaths)
    
    for oo in opaths:
        print(f'Fetching {oo}.')
    
    orealzs    = [Table.read(opath) for opath in opaths]
    mainreal   = orealzs[0]

    if len(realzs) == 1:
        print('Only one realization, nothing to be done. Exiting.')
        return mainreal

    keys = list(mainreal.meta.keys())

    for key in keys:
        print(key, mainreal.meta[key])

    if ('COLLATE' in keys) and (not force):
        print('Results have already been collated to zeroth realization. Exiting')
        return mainreal
    
    for col in ['FILLFACTOR', 'RAND_N8']:
        print(f'Solving for {col}')

        cols                     = [orealz[col].data.tolist() for orealz in orealzs]
        cols                     = np.array(cols).T
                
        mainreal[col]            = np.mean(cols, axis=1)
        
        # TODO: Check error on the mean.
        mainreal[col + '_STD']   =  np.std(cols, axis=1) / np.sqrt(len(realzs))
        
    print('Effective rand. density of {:.6f} to {:.6f}.'.format(mainreal.meta['RAND_DENS'], mainreal.meta['RAND_DENS'] * len(realzs)))

    mainreal.meta['COLLATE']     = 'TRUE'
    mainreal.meta['NREALZ']      = len(realzs)       

    fpath    = findfile(ftype='randoms', dryrun=dryrun, field=field, survey=survey, prefix=prefix, realz=0)
    boundary = Table.read(fpath, 'BOUNDARY')

    header   = fits.Header()

    hx       = fits.HDUList()
    hx.append(fits.PrimaryHDU(header=header))
    hx.append(fits.convenience.table_to_hdu(mainreal))
    hx.append(fits.convenience.table_to_hdu(boundary))
        
    if write:
        opath    = findfile(ftype='randoms_n8', dryrun=dryrun, field=field, survey=survey, prefix=prefix, realz=0)

        print(f'Writing {opath}.')

        hx.writeto(opath, overwrite=True)

    print('Finished.')
    
    return  mainreal

def process_one(run, pid=0, start=0.0):
    try:
        pid  = os.getpid()

    except Exception as e:
        print(e)

    bigtree  = run[0]
    comp     = run[1]

    msg      = 'POOL {}:  Creating {} tree for complement'.format(pid, len(comp))
    runtime  = calc_runtime(start, msg)

    # leafsize=5                                                                                                                                                                                        
    kd_tree  = KDTree(comp)

    msg      = 'POOL {}:  Querying {} tree for complement'.format(pid, len(bigtree.data))
    runtime  = calc_runtime(start, msg)

    # https://docs.scipy.org/doc/scipy/reference/generated/scipy.spatial.KDTree.query_ball_tree.html#scipy.spatial.KDTree.query_ball_tree                                                               
    # https://docs.scipy.org/doc/scipy/reference/generated/scipy.spatial.KDTree.count_neighbors.html#scipy.spatial.KDTree.count_neighbors                                                               
    indexes  = bigtree.query_ball_tree(kd_tree, r=sphere_radius)

    del kd_tree
    del bigtree
    del comp

    flat     = [len(idx) for idx in indexes]

    time.sleep(.1)

    return  flat

def fillfactor(log, field, dryrun, prefix, survey, oversample, nproc, realz, nooverwrite, debug=False):
    opath    = findfile(ftype='randoms_n8', dryrun=dryrun, field=field, survey=survey, prefix=prefix, realz=realz)

    if nooverwrite:
        overwrite_check(opath)

    start    = time.time()

    call_signature(dryrun, sys.argv)
    
    fpath          = findfile(ftype='randoms', dryrun=dryrun, field=field, survey=survey, prefix=prefix, oversample=oversample, realz=realz)
    overpoints_hdr = fitsio.read_header(fpath, ext=1)

    _overpoints    = fitsio.read(fpath, ext=1, columns=['CARTESIAN_X', 'CARTESIAN_Y', 'CARTESIAN_Z'])
    overpoints     = np.c_[_overpoints['CARTESIAN_X'], _overpoints['CARTESIAN_Y'], _overpoints['CARTESIAN_Z']]

    print(f'Fetching {fpath}')
    print('Fetched x{} oversampled randoms.'.format(overpoints_hdr['OVERSAMPLE']))

    # Read randoms file, split by field (DDP1, or not).                                                                                                                                                  
    # Note: realz handles oversampled realizations only.
    fpath          = findfile(ftype='randoms', dryrun=dryrun, field=field, survey=survey, prefix=prefix, realz=0)
    points_hdr     = fitsio.read_header(fpath, ext=1)

    _points        = fitsio.read(fpath, ext=1, columns=['CARTESIAN_X', 'CARTESIAN_Y', 'CARTESIAN_Z'])
    points         = np.c_[_points['CARTESIAN_X'], _points['CARTESIAN_Y'], _points['CARTESIAN_Z']]
 
    print(f'Fetching {fpath}')
    print('Fetched randoms of density {}.'.format(points_hdr['RAND_DENS']))
    
    del _points
    del _overpoints

    runtime        = calc_runtime(start, 'Reading {:.2f}M randoms'.format(len(overpoints) / 1.e6), xx=overpoints)

    print('Sorting.')

    idx            = np.argsort(points[:,0])
    points         = points[idx]

    idx            = np.argsort(overpoints[:,0])
    overpoints     = overpoints[idx]

    if debug:
        print('Assuming debug configuration.')

        debug_downsample = 5

        points           = points[::debug_downsample] 
        overpoints       = overpoints[::debug_downsample]

        print(f'Downsampling randoms by x{debug_downsample}.')

    runtime     = calc_runtime(start, 'Sorted randoms by X')

    # Chunked in x.
    split_idx   = np.arange(len(points))
    split_idx   = np.array_split(split_idx, 4 * nproc)

    nchunk      = len(split_idx)

    runs        = []

    for i, idx in enumerate(split_idx):
        split      = points[idx]

        xmin       = split[:,0].min()
        xmax       = split[:,0].max()
    
        buff       = .1  # [Mpc/h] 
    
        # Complement uses the oversampled version
        complement = (overpoints[:,0] > (xmin - sphere_radius - buff)) & (overpoints[:,0] < (xmax + sphere_radius + buff))
        complement = overpoints[complement]
    
        cmin       = complement[:,0].min()
        cmax       = complement[:,0].max() 

        print('{:d}\t{:.4f}\t{:.4f}\t{:.4f}\t{:.4f}\t{:d}\t{:d}'.format(i, xmin, xmax, cmin, cmax, len(split), len(complement)))

        # leafsize=5
        split_tree = KDTree(split)

        runs.append([split_tree, complement])

    runtime = calc_runtime(start, 'Created {} big trees and complement chunked by ...'.format(nchunk))

    del points
    del overpoints
    del split
    del split_idx

    runtime     = calc_runtime(start, 'Deleted rand.')

    runtime     = calc_runtime(start, 'POOL:  Counting < 8 Mpc/h pairs for small trees.')

    pool_start  = time.time()

    results     = [process_one(runs[0], pid=0, start=start)]

    split_time  = time.time() - pool_start
    split_time /= 60.

    runtime     = calc_runtime(start, 'POOL:  Expected runtime of {:.6f} minutes with {:d} proc. and split time {:.6f} mins'.format(nchunk * split_time / nproc, nproc, split_time))

    done_nsplit = 1

    # maxtasksperchild:  restart process after max tasks to contain resource leaks;
    with multiprocessing.get_context('spawn').Pool(nproc, maxtasksperchild=4) as pool:
        total   = (nchunk-1)

        with tqdm.tqdm(total=total) as pbar:
            for result in pool.imap(partial(process_one, start=start), iterable=runs[1:], chunksize=4):
                results.append(result)

                pbar.update()

                done_nsplit  += 1
                
                '''
                if (done_nsplit % nproc) == 0:
                    pool_time = (time.time() - pool_start)  / 60.
                    runtime   = calc_runtime(start, 'POOL:  New expected runtime of {:.3f} minutes with {:d} proc.'.format(nchunk * pool_time / done_nsplit, nproc))
                '''

        pool.close()

        # https://stackoverflow.com/questions/38271547/when-should-we-call-multiprocessing-pool-join                                                                                                       
        pool.join()

    runtime     = calc_runtime(start, 'POOL:  Done with queries of {} splits'.format(done_nsplit))
    # runtime   = calc_runtime(start, 'POOL:  Done with queries of {} splits with effective split time {}'.format(done_nsplit, pool_time / done_nsplit))

    flat_result = []
    
    for rr in results:
        flat_result += rr
 
    runtime                = calc_runtime(start, 'Reading randoms')

    fpath                  = findfile(ftype='randoms', dryrun=dryrun, field=field, survey=survey, prefix=prefix)
    rand                   = Table.read(fpath)

    runtime                = calc_runtime(start, 'Sorting randoms by X')
    
    rand.sort('CARTESIAN_X')

    if debug:
        rand               = rand[::debug_downsample] 

    # print(len(rand), len(flat_result))

    runtime                = calc_runtime(start, 'Assigning counts to randoms')

    rand['RAND_N8']        = np.array(flat_result).astype(np.int32)
    rand['FILLFACTOR']     = rand['RAND_N8'] / overpoints_hdr['NRAND8']

    rand.meta['RSPHERE']   = sphere_radius
    rand.meta['IMMUTABLE'] = 'FALSE'

    runtime                = calc_runtime(start, f'Reading {fpath}')

    boundary               = Table.read(fpath, 'BOUNDARY')
    
    header                 = fits.Header()

    hx                     = fits.HDUList()
    hx.append(fits.PrimaryHDU(header=header))
    hx.append(fits.convenience.table_to_hdu(rand))
    hx.append(fits.convenience.table_to_hdu(boundary))

    runtime                = calc_runtime(start, 'Writing {}.'.format(opath), xx=rand)

    hx.writeto(opath, overwrite=True)

    runtime                = calc_runtime(start, 'Finished')

    if log:
        sys.stdout.close()

    return 0


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Calculate fill factor using randoms.')
    parser.add_argument('--log', help='Create a log file of stdout.', action='store_true')
    parser.add_argument('-f', '--field', type=str, help='Select equatorial GAMA field: G9, G12, G15', default='G9')
    parser.add_argument('-d', '--dryrun', help='Dryrun.', action='store_true')
    parser.add_argument('-s', '--survey', help='Select survey.', default='gama')
    parser.add_argument('--prefix', help='filename prefix', default='randoms')
    parser.add_argument('--nproc', help='nproc', default=12, type=int)
    parser.add_argument('--realz', help='Realization number', default=0, type=np.int32)
    parser.add_argument('--nooverwrite',  help='Do not overwrite outputs if on disk', action='store_true')
    parser.add_argument('--oversample',   help='Oversampling factor for fillfactor counting.', default=2, type=int)
    parser.add_argument('--config',       help='Path to configuration file', type=str, default=findfile('config'))
    parser.add_argument('--debug', help='Trigger debug options, e.g. undersample', action='store_true')

    args        = parser.parse_args()
    log         = args.log
    field       = args.field.upper()
    dryrun      = args.dryrun
    prefix      = args.prefix
    survey      = args.survey.lower()
    oversample  = args.oversample
    debug       = args.debug

    # https://www.dur.ac.uk/icc/cosma/cosma5/                                                                                                                                                            
    nproc       = args.nproc
    realz       = args.realz
    nooverwrite = args.nooverwrite

    fields      = fetch_fields(survey)

    if log:
        logfile = findfile(ftype='randoms_n8', dryrun=False, field=field, survey=survey, prefix=prefix, realz=realz, log=True)

        print(f'Logging to {logfile}')

        sys.stdout = open(logfile, 'w')
    '''
    config = Configuration(args.config)
    config.update_attributes('fillfactor', args)
    config.write()
    '''
    assert field in fields, 'Error: Field not in fields'

    fillfactor(log, field, dryrun, prefix, survey, oversample, nproc, realz, nooverwrite, debug)
