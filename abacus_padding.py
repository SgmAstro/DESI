import argparse
from astropy.table import vstack, Table

# HACK to deal with path issue
import os
import sys
sys.path.append('{}/DESI'.format(os.environ['HOME']))

from findfile      import findfile, fetch_fields, overwrite_check, write_desitable


def edge_padding(mock, opath, dryrun, pad=10):
    '''
    Inputs the abacus mock and adds the appropriate galaxies to remove the edge.
    (Specifically we reflect on the edges).
    
    # TODO: work out a more computationally efficient way of doing this.
    
    # e.g: work in mod 360 for left and right sides.
    '''
    boundary = 1000 # upper bound on x, y, z
    
    if dryrun:
        pad = 5
        boundary = 20
    
    if not dryrun:
        min_x, max_x = min(mock['CARTESIAN_X']), max(mock['CARTESIAN_X'])
        min_y, max_y = min(mock['CARTESIAN_Y']), max(mock['CARTESIAN_Y'])
        min_z, max_z = min(mock['CARTESIAN_Z']), max(mock['CARTESIAN_Z'])

        assert (min_x < -999) and (max_x > 999), f'Not a cubic box, edge correction may fail here.'
        assert (min_y < -999) and (max_y > 999), f'Not a cubic box, edge correction may fail here.'
        assert (min_z < -999) and (max_z > 999), f'Not a cubic box, edge correction may fail here.'

    # mark all galaxies in mock as 'real' for LF purposes
    mock['REP_GAL'] = 0
    
    print('Padding left and right sides.')
        
    # add left and right padding
    
    # add left and right padding
    r_side = mock[mock['CARTESIAN_X'] > boundary-pad]
    l_side = mock[mock['CARTESIAN_X'] < pad-boundary]

    r_side['CARTESIAN_X'] = r_side['CARTESIAN_X'] + 2*(boundary - r_side['CARTESIAN_X'])
    l_side['CARTESIAN_X'] = l_side['CARTESIAN_X'] - 2*(boundary + l_side['CARTESIAN_X'])

    l_side['REP_GAL'] = 1
    r_side['REP_GAL'] = 1

    mock = vstack([mock, l_side, r_side])
    
    print('Padding top and bottom sides.')

    
    # add top and bottom padding
    u_side = mock[mock['CARTESIAN_Y'] > boundary-pad]
    d_side = mock[mock['CARTESIAN_Y'] < pad-boundary]

    u_side['CARTESIAN_Y'] = u_side['CARTESIAN_Y'] + 2*(boundary - u_side['CARTESIAN_Y'])
    d_side['CARTESIAN_Y'] = d_side['CARTESIAN_Y'] - 2*(boundary + d_side['CARTESIAN_Y'])
    
    u_side['REP_GAL'] = 1
    d_side['REP_GAL'] = 1

    mock = vstack([mock, u_side, d_side])

    print('Padding front and back sides.')

    
    # add front and back padding
    f_side = mock[mock['CARTESIAN_Z'] > boundary-pad]
    b_side = mock[mock['CARTESIAN_Z'] < pad-boundary]

    f_side['CARTESIAN_Z'] = f_side['CARTESIAN_Z'] + 2*(boundary - f_side['CARTESIAN_Z'])
    b_side['CARTESIAN_Z'] = b_side['CARTESIAN_Z'] - 2*(boundary + b_side['CARTESIAN_Z'])
    
    f_side['REP_GAL'] = 1
    b_side['REP_GAL'] = 1

    mock = vstack([mock, f_side, b_side])
    
    print('Writing {}.'.format(opath))
    
    write_desitable(opath, mock)
    
    

if __name__ == '__main__':

    parser  = argparse.ArgumentParser(description='Gen abacus gold.')
    parser.add_argument('-d', '--dryrun', help='Dryrun.', action='store_true')
    parser.add_argument('--nooverwrite',  help='Do not overwrite outputs if on disk', action='store_true')
            
    args        = parser.parse_args()
    dryrun      = args.dryrun
    nooverwrite = args.nooverwrite
    
    
    # TODO: change ftype to a later file    
    fpath = findfile(ftype='gold', survey='abacus', dryrun=dryrun)
    
    opath = fpath.replace('gold', 'gold_padded')
    
    print('Reading in mock.')

    mock = Table.read(fpath)
    
    edge_padding(mock, opath, dryrun=dryrun, pad=10)