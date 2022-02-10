import numpy         as     np

from   astropy.table import Table, unique
from   ros_tools     import tile2rosette, calc_rosr


def desi_randoms(field):
    # TODO: Handle NERSC (generation) vs Cosma (read only).
    # Randoms uniform on the sphere with density 2500 per sq. deg., available to an assigned fiber.  
    rand = Table.read('/global/cfs/cdirs/desi/survey/catalogs/SV3/LSS/random0/rancomb_brightwdup_Alltiles.fits')
    # rand.pprint()

    # TODO:  Check TARGETID is a unique idetifier, or bug.  If not, use RA. 
    rand             = unique(rand, keys='TARGETID')
    rand['ROS']      = tile2rosette(rand['TILEID'])
    rand['ROS_DIST'] = 1.e99

    rand = rand[rand['ROS'] == int(field[1:])]
    
    for rosn in np.unique(rand['ROS']):
        isin = (rand['ROS'].data == rosn)

        new_dist = calc_rosr(rosn, rand['RA'][isin], rand['DEC'][isin])
    
        rand['ROS_DIST'][isin] = np.minimum(rand['ROS_DIST'][isin], new_dist)

    # rand.pprint()

    # TODO:  check.
    rand = rand[(rand['ROS_DIST'].data > 0.5) & (rand['ROS_DIST'].data < 1.5)]
    
    rand.rename_column('RA',  'RANDOM_RA')
    rand.rename_column('DEC', 'RANDOM_DEC')
    
    return  rand
