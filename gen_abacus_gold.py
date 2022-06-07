import numpy as np
import matplotlib.pyplot as plt
import h5py

from astropy.table import Table


# HACK to deal with path issue
import os
import sys
sys.path.append('{}/DESI'.format(os.environ['HOME']))


from findfile      import findfile, fetch_fields, overwrite_check, write_desitable


def abacus_gold():
    
    path = "/project/projectdirs/desi/cosmosim/FirstGenMocks/AbacusSummit/CutSky/BGS/z0.200/"
    opath = findfile(ftype='gold', survey='abacus')

    f = h5py.File(path+"cutsky_BGS_z0.200_AbacusSummit_base_c000_ph000.hdf5", "r")

    # TODO: clean up with f['Data'].keys()
    stat    = f['Data/STATUS'][...] 
    Mr      = f["Data/abs_mag"][...]
    r       = f["Data/app_mag"][...]
    ra      = f["Data/ra"][...]
    dec     = f["Data/dec"][...]
    gmr     = f["Data/g_r"][...]
    gmr_obs = f["Data/g_r_obs"][...]
    z_obs   = f["Data/z_obs"][...]
    z_cos   = f["Data/z_cos"][...]
    galtype = f["Data/galaxy_type"][...]
    halo    = f["Data/halo_mass"][...]

    f.close()

    mock  = Table(np.c_[ra, dec, gmr, gmr_obs, z_obs, r, Mr, stat, galtype, halo], names=['RA', 'DEC', 'G-R', 'G-R_OBS', 'Z_OBS', 'R_MAG', 'M', 'STATUS', 'GALTYPE', 'M_HALO'])

    print('Writing {}.'.format(opath))
    
    write_desitable(opath, mock)


if __name__ == '__main__':

    abacus_gold()