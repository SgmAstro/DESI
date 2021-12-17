import os
import fitsio
import numpy as np
import matplotlib.pyplot as plt

from   astropy.table import Table
from   scipy.spatial import KDTree
from   cartesian import cartesian


fpath = os.environ['CSCRATCH'] + '/norberg/GAMA4/gama_gold_ddp.fits'

dat = fitsio.read(fpath)
dat = Table(dat)

# dat = dat[:5000]

# add cartesian coordinates to catalogue
xyz = cartesian(dat['RA'], dat['DEC'], dat['ZGAMA'])

dat['CARTESIAN_X'] = xyz[:,0]
dat['CARTESIAN_Y'] = xyz[:,1]
dat['CARTESIAN_Z'] = xyz[:,2]

points      = np.c_[dat['CARTESIAN_X'], dat['CARTESIAN_Y'], dat['CARTESIAN_Z']]
points      = np.array(points, copy=True)

kd_tree_all  = KDTree(points)

# randoms.
field  = 'G9' 
realz  = 0

rpath  = os.environ['CSCRATCH'] + '/desi/BGS/Sam/randoms_bd_{}_{:d}.fits'.format(field, realz)
rand   = fitsio.read(rpath)

# rand   = rand[:100000]

rpoints = np.c_[rand['CARTESIAN_X'], rand['CARTESIAN_Y'], rand['CARTESIAN_Z']]
rpoints = np.array(rpoints, copy=True)

print('Creating big rand. tree.')

big_tree = KDTree(rpoints)

print('Querying tree for closest rand.')

dd, ii   = big_tree.query([x for x in points], k=1)

dat['RANDSEP'] = dd
dat['RANDMATCH'] = rand['RANDID'][ii]
dat['BOUND_DIST'] = rand['BOUND_DIST'][ii]

# TODO: Get from header.
nrand8 = 1072.3302924

dat['FILLFACTOR'] = rand['N8'][ii] / nrand8

for idx in range(3):
    ddp_idx      = idx + 1
    
    ddp          = dat[dat['DDP'][:,idx] == 1]
    points_ddp   = np.c_[ddp['CARTESIAN_X'], ddp['CARTESIAN_Y'], ddp['CARTESIAN_Z']]
    points_ddp   = np.array(points_ddp, copy=True)

    print('Building tree for DDP {}'.format(ddp_idx))
    
    kd_tree_ddp  = KDTree(points_ddp)

    print('Querying tree for DDP {}'.format(ddp_idx))

    indexes_ddp  = kd_tree_all.query_ball_tree(kd_tree_ddp, r=8.)

    dat['DDP{:d}_N8'.format(ddp_idx)] = np.array([len(idx) for idx in indexes_ddp])

dat.pprint()

# TODO: Inherit the DDP header and propagate. 
meta     = {'DDP1_ZMIN': 0.039069999009370804,\
            'DDP1_ZMAX': 0.2483299970626831,\
            'DDP1_VZ': 6451530.309761727,\
            'DDP1_DENS': 0.005383528919866882,\
            'DDP2_ZMIN': 0.03914999961853027,\
            'DDP2_ZMAX': 0.18308000266551971,\
            'DDP2_VZ': 2679079.7557868413,\
            'DDP2_DENS': 0.009928035902084674,\
            'DDP3_ZMIN': 0.03903000056743622,\
            'DDP3_ZMAX': 0.09973999857902527,\
            'DDP3_VZ': 432372.2344703941,\
            'DDP3_DENS': 0.018396185892331243,\
            'VOL8': (4./3.)*np.pi*(8.**3.)}

for x in meta.keys():
    print('{}\t\t{}'.format(x.ljust(20), meta[x]))

print('Writing {}'.format(fpath.replace('ddp', 'ddp_n8')))

dat.write(fpath.replace('ddp', 'ddp_n8'), overwrite=True, format='fits')
