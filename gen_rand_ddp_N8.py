import os
import fitsio
import numpy as np
import matplotlib.pyplot as plt
import argparse

from astropy.table import Table
from scipy.spatial import KDTree
from cartesian import cartesian
from delta8_limits import dd8_limits, delta8_tier


parser = argparse.ArgumentParser(description='Calculate DDP1 N8 for all randoms.')
parser.add_argument('-f', '--field', type=str, help='Select equatorial GAMA field: G9, G12, G15', required=True)
parser.add_argument('-d', '--dryrun', help='Dryrun.', action='store_true')

args   = parser.parse_args()
field  = args.field.upper()
dryrun = args.dryrun

realz  = 0

fpath  = os.environ['GOLD_DIR'] + '/gama_gold_ddp.fits'

if dryrun:
    fpath = fpath.replace('.fits', '_dryrun.fits')

dat = Table.read(fpath)

fpath = os.environ['RANDOMS_DIR'] + '/randoms_bd_{}_{:d}.fits'.format(field, realz)

if dryrun:
    fpath = fpath.replace('.fits', '_dryrun.fits')

rand = Table.read(fpath)

# Propagate header 'DDP1_ZMIN' etc. to randoms.
rand.meta.update(dat.meta)

points       = np.c_[rand['CARTESIAN_X'], rand['CARTESIAN_Y'], rand['CARTESIAN_Z']]
points       = np.array(points, copy=True)

kd_tree_rand = KDTree(points)

for idx in range(3):
    ddp_idx      = idx + 1

    print('Solving for DDP {}'.format(ddp_idx))
    
    ddp          = dat[dat['DDP'][:,idx] == 1]
    points_ddp   = np.c_[ddp['CARTESIAN_X'], ddp['CARTESIAN_Y'], ddp['CARTESIAN_Z']]
    points_ddp   = np.array(points_ddp, copy=True)
    
    kd_tree_ddp  = KDTree(points_ddp)
    
    indexes_ddp  = kd_tree_rand.query_ball_tree(kd_tree_ddp, r=8.)

    rand['DDP{:d}_N8'.format(ddp_idx)] = np.array([len(idx) for idx in indexes_ddp])


rand['DDP1_DELTA8'] = (rand['DDP1_N8'] / (dat.meta['VOL8'] * dat.meta['DDP1_DENS']) / rand['FILLFACTOR']) - 1.
rand['DDP2_DELTA8'] = (rand['DDP2_N8'] / (dat.meta['VOL8'] * dat.meta['DDP2_DENS']) / rand['FILLFACTOR']) - 1.
rand['DDP3_DELTA8'] = (rand['DDP3_N8'] / (dat.meta['VOL8'] * dat.meta['DDP3_DENS']) / rand['FILLFACTOR']) - 1.

rand['DDP1_DELTA8_TIER'] = delta8_tier(rand['DDP1_DELTA8'])

utiers = np.unique(rand['DDP1_DELTA8_TIER'].data)

ddp1_zmin = dat.meta['DDP1_ZMIN']
ddp1_zmax = dat.meta['DDP1_ZMAX']

print('Unique tiers: {}'.format(utiers))

print('Found redshift limits: {:.3f} < z < {:.3f}'.format(ddp1_zmin, ddp1_zmax))

for ut in utiers:    
    ddp1_rand = rand[(rand['Z'] > ddp1_zmin) & (rand['Z'] < ddp1_zmax)]    
    in_tier   = (ddp1_rand['DDP1_DELTA8_TIER'].data == ut)
        
    rand.meta['DDP1_d{}_VOLFRAC'.format(ut)]   = np.mean(in_tier)
    rand.meta['DDP1_d{}_TIERMEDd8'.format(ut)] = np.median(ddp1_rand['DDP1_DELTA8'].data[in_tier])

    print('DDP1_d{}_VOLFRAC OF {:.4f} ADDED.'.format(ut, np.mean(in_tier)))
    print('DDP1_d{}_TIERMEDd8 OF {:.4f} ADDED.'.format(ut, rand.meta['DDP1_d{}_TIERMEDd8'.format(ut)]))
    
opath = fpath.replace('bd', 'bd_ddp_n8')
    
print('Writing {}'.format(opath))

rand.write(opath, format='fits', overwrite=True)
