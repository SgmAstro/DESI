import os
import sys
import fitsio
import argparse
import runtime
import numpy as np
import matplotlib.pyplot as plt

from   astropy.table import Table, vstack
from   scipy.spatial import KDTree
from   delta8_limits import delta8_tier, d8_limits
from   findfile      import findfile, fetch_fields, overwrite_check, gather_cat, write_desitable, fetch_header
from   config        import Configuration
from   bitmask       import lumfn_mask, consv_mask, update_bit
from   delta8_limits import d8_limits
from   runtime       import calc_runtime
from   params        import fillfactor_threshold, oversample_nrealisations, sphere_radius

parser = argparse.ArgumentParser(description='Generate DDP1 N8 for all gold galaxies.')
parser.add_argument('--log', help='Create a log file of stdout.', action='store_true')
parser.add_argument('-d', '--dryrun', help='Dryrun.', action='store_true')
parser.add_argument('-s', '--survey', help='Select survey', default='gama')
parser.add_argument('--realz', help='Realization', default=0, type=int)
parser.add_argument('--oversample', help='Oversample', default=2, type=int)
parser.add_argument('--oversample_nrealisations', help='Oversample realization number', default=None)
parser.add_argument('--nooverwrite',  help='Do not overwrite outputs if on disk', action='store_true')

args        = parser.parse_args()
log         = args.log
realz       = args.realz
dryrun      = args.dryrun
survey      = args.survey.lower()
oversample  = args.oversample

if args.oversample_nrealisations != None:
    oversample_nrealisations = int(args.oversample_nrealisations)
  
    print(f'Overriding number of oversampled realizations used with {oversample_nrealisations}')

fields      = fetch_fields(survey)

fpath       = findfile(ftype='ddp',    dryrun=dryrun, survey=survey)
opath       = findfile(ftype='ddp_n8', dryrun=dryrun, survey=survey)

if log:
    logfile = findfile(ftype='ddp_n8', dryrun=False, survey=survey, log=True)

    print(f'Logging to {logfile}')
        
    sys.stdout = open(logfile, 'w')

if args.nooverwrite:
    overwrite_check(opath)
    
# Read ddp cat.    
dat           = Table.read(fpath)

print('Reading: {} with length {}'.format(fpath, len(dat)))

assert 'DDP1_DENS' in dat.meta

points       = np.c_[dat['CARTESIAN_X'], dat['CARTESIAN_Y'], dat['CARTESIAN_Z']]
points       = np.array(points, copy=True)

kd_tree_all  = KDTree(points)

# Oversampled randoms 
prefix           = 'randoms_ddp1'
dat['RAND_N8']   = 0.

for realz in np.arange(oversample_nrealisations):
    print(f'\n\nSolving for galaxy fillfactors with oversampled realization {realz}.')

    rpaths       = [findfile(ftype='randoms', dryrun=dryrun, field=ff, survey=survey, prefix=prefix, oversample=oversample, realz=realz) for ff in fields]

    for rpath in rpaths:
        print('Fetching: {}'.format(rpath))

    orand        = gather_cat(rpaths)

    orpoints     = np.c_[orand['CARTESIAN_X'], orand['CARTESIAN_Y'], orand['CARTESIAN_Z']]

    print('Creating oversample rand. tree.')

    obig_tree       = KDTree(orpoints)
    
    indexes_dat     = kd_tree_all.query_ball_tree(obig_tree, r=8.)
    dat['RAND_N8'] += np.array([len(idx) for idx in indexes_dat])

    print('After solving for realization {}, median number of randoms per 8-sphere is {}'.format(realz, np.median(dat['RAND_N8'])))
    
del orand
del orpoints
del obig_tree

hpath               = findfile(ftype='randoms_n8', dryrun=dryrun, field=fields[0], survey=survey, prefix=prefix, oversample=1, realz=0)

print(f'Fetching header information from {hpath}')

onrand8             = oversample_nrealisations * oversample * fetch_header(fpath=hpath, name='NRAND8')
ordens              = oversample_nrealisations * oversample * fetch_header(fpath=hpath, name='RAND_DENS') 

dat['FILLFACTOR']   = dat['RAND_N8'] / onrand8

print('Normalised galaxy fill factors with {:.2f} expected randoms per 8-sphere (density: {:.6e}).'.format(onrand8, ordens))


# ----  Find closest matching oversampled random to inherit bounddist  ----
print('Finding bound dist measure.')

bpaths              = [findfile(ftype='randoms_n8', dryrun=dryrun, field=ff, survey=survey, prefix=prefix) for ff in fields]
boundary            = [Table.read(bpath, 'BOUNDARY') for bpath in bpaths]

# TODO Note: BOUNDID will not be unique.
boundary            = vstack(boundary)
boundary            = np.c_[boundary['CARTESIAN_X'], boundary['CARTESIAN_Y'], boundary['CARTESIAN_Z']]
boundary_tree       = KDTree(boundary)

body                = np.c_[dat['CARTESIAN_X'], dat['CARTESIAN_Y'], dat['CARTESIAN_Z']]
split               = [x for x in body]

dd, ii              = boundary_tree.query(split, k=1)
dat['BOUND_DIST']   = dd

dat['FILLFACTOR'][dat['BOUND_DIST'] > sphere_radius] = 1.


# ----  Find closest matching random to inherit fill factor  ----
# Read randoms bound_dist.
rpaths              = [findfile(ftype='randoms_bd', dryrun=dryrun, field=ff, survey=survey, prefix=prefix, oversample=1, realz=0) for ff in fields]

for rpath in rpaths:
    print('Reading: {}'.format(rpath))

rand                = gather_cat(rpaths)

print('Retrieved galaxies for {}'.format(np.unique(dat['FIELD'].data)))
print('Retrieved randoms for {}'.format(np.unique(rand['FIELD'].data)))

for i, rpath in enumerate(rpaths):
    dat.meta['RPATH_{}'.format(i)] = rpath

rpoints  = np.c_[rand['CARTESIAN_X'], rand['CARTESIAN_Y'], rand['CARTESIAN_Z']]

print('Creating big rand. tree.')

big_tree = KDTree(rpoints)

print('Querying tree for closest rand.')

dd, ii   = big_tree.query([x for x in points], k=1)

# Find closest random for bound_dist and fill factor. 
# These randoms are split by field.
dat['rRANDSEP']    = dd
dat['rRANDMATCH']  = rand['RANDID'][ii]
dat['rBOUND_DIST'] = rand['BOUND_DIST'][ii]
dat['rFILLFACTOR'] = rand['FILLFACTOR'][ii]

update_bit(dat['IN_D8LUMFN'], lumfn_mask, 'FILLFACTOR', dat['FILLFACTOR'].data < fillfactor_threshold)

if not dryrun:
    match_sep = 6.5

    # Typically, bounded by 1.6
    # assert  np.all(dat['rRANDSEP'].data < match_sep), 'Failed to find matching random with < 5 Mpc/h separation.'

    if not np.all(dat['rRANDSEP'].data < match_sep):
        # Note: DESI randoms are less dense, larger expected separation.
        print('WARNING: poor random match, with maximum comoving random separation >3Mpc/h.')

        poor_match = dat['rRANDSEP'].data > match_sep

        print(dat[poor_match])

# ----  Calculate DDPX_N8 for each gama gold galaxy.  ----
for idx in range(3):
    # Calculate DDP1/2/3 N8 for all gold galaxies.
    ddp_idx      = idx + 1

    dat['DDP{:d}_N8'.format(ddp_idx)] = -99
    
    for field in fields:
        print('Building tree for DDP {} and field {}'.format(ddp_idx, field))

        in_field      = dat['FIELD'] == field
        dat_field     = dat[in_field]

        ddp           = dat_field[dat_field['DDP'][:,idx] == 1]
        points_ddp    = np.c_[ddp['CARTESIAN_X'], ddp['CARTESIAN_Y'], ddp['CARTESIAN_Z']]
        points_ddp    = np.array(points_ddp, copy=True)
        
        kd_tree_ddp   = KDTree(points_ddp)

        print('Querying tree for DDP {}'.format(ddp_idx))

        indexes_ddp   = kd_tree_all.query_ball_tree(kd_tree_ddp, r=8.)

        counts        = np.array([len(idx) for idx in indexes_ddp]) 

        dat['DDP{:d}_N8'.format(ddp_idx)][in_field] = counts[in_field] 

##  Derived.
dat.meta['VOL8']   = (4./3.)*np.pi*(8.**3.)

dat['DDP1_DELTA8'] = ((dat['DDP1_N8'] / (dat.meta['VOL8'] * dat.meta['DDP1_DENS']) / dat['FILLFACTOR'])) - 1. 

##  
outwith = (dat['ZSURV'] > dat.meta['DDP1_ZMIN']) & (dat['ZSURV'] < dat.meta['DDP1_ZMAX'])
outwith = ~outwith

if not dryrun:
    # Insufficient randoms in a dryrun.
    outwith = outwith | (dat['FILLFACTOR']  < fillfactor_threshold)

dat['DDP1_DELTA8'][outwith] = -99.
dat['DDP1_DELTA8_TIER']     = delta8_tier(dat['DDP1_DELTA8'])

dat.pprint()

# TODO: Check
if 'ddp1' not in prefix:
    dat['DDP2_DELTA8'] = ((dat['DDP2_N8'] / (dat.meta['VOL8'] * dat.meta['DDP2_DENS']) / dat['FILLFACTOR'])) - 1. 
    dat['DDP3_DELTA8'] = ((dat['DDP3_N8'] / (dat.meta['VOL8'] * dat.meta['DDP3_DENS']) / dat['FILLFACTOR'])) - 1. 

for x in dat.meta.keys():
    print('{}\t\t{}'.format(x.ljust(20), dat.meta[x]))

print('Writing {}'.format(opath))

write_desitable(opath, dat)

#  ----  Generate ddp_n8_d0 files for LF(d8) files, limited to DDP1 (and redshift range)  ----
dat                     = dat[(dat['ZSURV'] > dat.meta['DDP1_ZMIN']) & (dat['ZSURV'] < dat.meta['DDP1_ZMAX'])]
dat['DDP1_DELTA8_TIER'] = delta8_tier(dat['DDP1_DELTA8'])

utiers                  = np.unique(dat['DDP1_DELTA8_TIER'].data)

if -99 in utiers:
    utiers = utiers.tolist()    
    utiers.remove(-99)
    utiers = np.array(utiers)

for ii, xx in enumerate(d8_limits):
    dat.meta['D8{}LIMS'.format(ii)] = str(xx)

if not np.all(np.isin(np.arange(9), utiers)):
    print('WARNING: MISSING d8 TIERS ({})'.format(utiers))
    
else:
    print(utiers)

print('Delta8 spans {:.4f} to {:.4f} over {} tiers.'.format(dat['DDP1_DELTA8'].min(), dat['DDP1_DELTA8'].max(), utiers))

for tier in np.arange(len(d8_limits)):
    print()
    print('---- d{} ----'.format(tier))

    isin     = (dat['DDP1_DELTA8_TIER'].data == tier)    
    to_write = dat[isin]

    dat.meta['DDP1_D{}_NGAL'.format(tier)] = len(to_write)

    assert 'AREA' in dat.meta.keys()
    assert 'AREA' in to_write.meta.keys()

    print('Available fields in tier: {}'.format(np.unique(dat['FIELD'].data)))

    for field in fields:    
        isin           = to_write['FIELD'] == field
        to_write_field = to_write[isin]

        opath_field    = findfile('ddp_n8_d0', dryrun=dryrun, field=field, utier=tier, survey=survey, realz=realz)  

        print('Writing {} galaxies from field {} to {}.'.format(len(to_write_field), np.unique(to_write_field['FIELD'].data), opath_field))

        to_write_field.meta['AREA'] = to_write.meta['AREA'] / len(fields)

        write_desitable(opath_field, to_write_field)

print('\n\nDone.\n\n')

if log:
    sys.stdout.close()
