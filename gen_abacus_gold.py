import numpy as np
import matplotlib.pyplot as plt
import h5py
import argparse

from astropy.table import Table


# HACK to deal with path issue
import os
import sys
sys.path.append('{}/DESI'.format(os.environ['HOME']))


from findfile      import findfile, fetch_fields, overwrite_check, write_desitable


def abacus_gold(dryrun=False):
    
    #path = "/project/projectdirs/desi/cosmosim/FirstGenMocks/AbacusSummit/CutSky/BGS/z0.200/"
    path = "/project/projectdirs/desi/cosmosim/FirstGenMocks/AbacusSummit/CubicBox/BGS/z0.200/AbacusSummit_base_c000_ph000/"

    opath = findfile(ftype='gold', survey='abacus', dryrun=dryrun)

    f = h5py.File(path+"BGS_box_ph000.hdf5", "r")
    
    
    # TODO: clean up with f['Data'].keys()
    pos = f["Data/pos"][...] # position, in comoving Mpc/h, in range -1000 < x < 1000 Mpc/h
    #vel = f["Data/vel"][...]
    abs_mag = f["Data/abs_mag"][...]
    g_r = f["Data/g_r"][...]
    #hmass = f["Data/halo_mass"][...]
    galtype = f["Data/galaxy_type"][...]
    f.close()

    abacus_gold = Table([pos[:,0], pos[:,1], pos[:,2], abs_mag, g_r, galtype], names=['CARTESIAN_X', 'CARTESIAN_Y', 'CARTESIAN_Z', 'DETMAG', 'GMR', 'GALTYPE'])
    
    print(len(abacus_gold))
    
    
    # restrict magnitudes
    '''
    abacus_gold.meta['RMAX'] = 12.0
    abacus_gold.meta['RLIM'] = 19.8

    isin        = (abacus_gold.meta['RMAX'] <= abacus_gold['DETMAG']) & (abacus_gold['DETMAG'] <= abacus_gold.meta['RLIM'])
    abacus_gold = abacus_gold[isin]
    
    print(len(abacus_gold))
    '''
        
    if dryrun:
        isin = (abacus_gold['CARTESIAN_X'] <= 20) & (abacus_gold['CARTESIAN_X'] >= -20) & (abacus_gold['CARTESIAN_Y'] <= 20) & (abacus_gold['CARTESIAN_Y'] >= -20) &  (abacus_gold['CARTESIAN_Z'] <= 20) & (abacus_gold['CARTESIAN_Z'] >= -20)
        
        abacus_gold = abacus_gold[isin]
    
    print(len(abacus_gold))

    
    print('Writing {}.'.format(opath))
    
    write_desitable(opath, abacus_gold)


if __name__ == '__main__':

    parser  = argparse.ArgumentParser(description='Gen abacus gold.')
    parser.add_argument('-d', '--dryrun', help='Dryrun.', action='store_true')
    parser.add_argument('--nooverwrite',  help='Do not overwrite outputs if on disk', action='store_true')
            
    args        = parser.parse_args()
    dryrun      = args.dryrun
    nooverwrite = args.nooverwrite
    
    abacus_gold(dryrun)