import numpy         as     np

from   astropy.table import Table, unique
from   ros_tools     import tile2rosette, calc_rosr


def desi_randoms():
    rand = Table.read('/global/cfs/cdirs/desi/survey/catalogs/SV3/LSS/random0/rancomb_brightwdup_Alltiles.fits')
    # rand.pprint()

    rand             = unique(rand, keys='TARGETID')
    rand['ROS']      = tile2rosette(rand['TILEID'])
    rand['ROS_DIST'] = 1.e99

    for rosn in np.unique(rand['ROS']):
        isin = (rand['ROS'].data == rosn)

        new_dist = calc_rosr(rosn, rand['RA'][isin], rand['DEC'][isin])
    
        rand['ROS_DIST'][isin] = np.minimum(rand['ROS_DIST'][isin], new_dist)

    # rand.pprint()

    rand.rename_column('RA',  'RANDOM_RA')
    rand.rename_column('DEC', 'RANDOM_DEC')
    
    return  rand
