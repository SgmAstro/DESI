from astropy.table import vstack

# HACK to deal with path issue
import os
import sys
sys.path.append('{}/DESI'.format(os.environ['HOME']))

from findfile      import findfile, fetch_fields, overwrite_check, write_desitable


def edge_padding(mock, opath, pad=10):
    '''
    Inputs the abacus mock and adds the appropriate galaxies to remove the edge.
    (Specifically we reflect on the edges).
    
    # TODO: work out a more computationally efficient way of doing this.
    
    # e.g: work in mod 360 for left and right sides.
    '''
    
    min_ra, max_ra = min(mock['RA']), max(mock['RA'])
    min_dec, max_dec = min(mock['DEC']), max(mock['DEC'])
    
    assert (min_ra < 1) and (max_ra > 359), f'Not an all-sky mock, edge correction will fail here.'
    assert (min_dec < -89) and (max_dec > 89), f'Not an all-sky mock, edge correction will fail here.'
    
    # mark all galaxies in mock as 'real' for LF purposes
    mock['REP_GAL'] = 0
    
    print('Padding left and right sides.')
        
    # add left and right padding
    l_side = mock[mock['RA'] > 360-pad]
    r_side = mock[mock['RA'] < pad]

    l_side['RA'] = l_side['RA'] - 360
    r_side['RA'] = r_side['RA'] + 360

    l_side['REP_GAL'] = 1
    r_side['REP_GAL'] = 1

    mock = vstack([mock, l_side, r_side])
    
    print('Padding top and bottom sides.')

    
    # add top and bottom padding
    u_side = mock[mock['DEC'] > 90-pad]
    d_side = mock[mock['DEC'] < pad-90]

    u_side['DEC'] = (90 - u_side['DEC']) + 90
    d_side['DEC'] = (-90 - d_side['DEC']) - 90

    # for convenience, work with positive integers
    mock['RA'] += pad
    u_side['RA'] += pad
    d_side['RA'] += pad
    
    u_side['RA'] = (u_side['RA'] + 180) % (360 + 2*pad)
    d_side['RA'] = (d_side['RA'] + 180) % (360 + 2*pad)

    u_side['REP_GAL'] = 1
    d_side['REP_GAL'] = 1

    mock = vstack([mock, u_side, d_side])
    
    mock['RA'] = mock['RA'] - pad
    
    print('Writing {}.'.format(opath))
    
    write_desitable(opath, mock)
    
    

if __name__ == '__main__':

    # TODO: add findfile
    fpath = '/global/cscratch1/sd/ldrm11/desi/abacus/abacus_gold.fits'
    opath = fpath.replace('gold.fits', 'gold_padded.fits')
    
    mock = Table.read(fpath)
    
    edge_padding(mock, opath, pad=10)