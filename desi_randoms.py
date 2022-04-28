import os
import numpy         as     np

from   bitmask       import lumfn_mask
from   astropy.table import Table, unique, vstack
from   ros_tools     import tile2rosette, calc_rosr


def desi_randoms(ros, nrealz=15, dryrun=False):
    assert  'NERSC_HOST' in os.environ.keys()

    # Randoms uniform on the sphere with density 2500 per sq. deg., available to an assigned fiber.      
    rpaths = ['/global/cfs/cdirs/desi/survey/catalogs/SV3/LSS/random{}/rancomb_brightwdup_Alltiles.fits'.format(idx) for idx in range(1 + nrealz)]
    rand   = vstack([Table.read(xx) for xx in rpaths])
                      
    # TODO:  Check TARGETID is a unique identifier, or bug.  If not, use RA. 
    rand             = unique(rand, keys='TARGETID')
    rand['ROS']      = tile2rosette(rand['TILEID'])

    rand['ROS_DIST'] = 1.e99

    rand             = rand[rand['ROS'] == ros]
    
    for rosn in np.unique(rand['ROS']):
        isin = (rand['ROS'].data == rosn)

        new_dist = calc_rosr(rosn, rand['RA'][isin], rand['DEC'][isin])
    
        rand['ROS_DIST'][isin] = np.minimum(rand['ROS_DIST'][isin], new_dist)

    # rand.pprint()

    # TODO:  check.
    if dryrun:
        limits = [0.9, 1.1]
    else:
        limits = [0.5, 1.5]

    hi_comp             = (rand['ROS_DIST'].data > limits[0]) & (rand['ROS_DIST'].data < limits[1])
    rand['IN_D8LUMFN']  = ~hi_comp * lumfn_mask.DESI_HICOMP
    
    rand.rename_column('RA',  'RANDOM_RA')
    rand.rename_column('DEC', 'RANDOM_DEC')

    if dryrun:
        rand = rand[rand['IN_D8LUMFN'] == 0]

    # Must come after dryrun.
    rand.meta['AREA'] = len(rand) / 2500. / nrealz
    rand.meta['NRAND'] = len(rand)
    rand.meta['IMMUTABLE'] = 'TRUE'
        
    return  rand
