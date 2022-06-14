import sys
import h5py
import numpy               as     np

from   cosmo               import cosmo
from   cartesian           import cartesian
from   astropy.table       import Table
from   ros_tools           import tile2rosette
from   rotate_mock         import rotate_mock


# const. luminosity function.
fpath   = '/cosma5/data/durham/tccb49/MXXL_Lightcone/Final_Mock/v5/fixed_luminosity_function/full_sky/galaxy_catalogue_z0p2.hdf5'

f       = h5py.File(fpath, mode='r')

ra      = f["Data/ra"][...]
dec     = f["Data/dec"][...]
zobs    = f["Data/z_obs"][...]
z       = f["Data/z_cos"][...]
r       = f["Data/app_mag"][...]
M       = f["Data/abs_mag"][...]
gmr     = f["Data/g_r"][...]
tt      = f["Data/galaxy_type"][...]
hmass   = f["Data/halo_mass"][...]
mxxlid  = f["Data/mxxl_id"][...]

isin    = (r <= 19.5)

print(len(ra), len(nmock), 100. * np.mean(isin))

exit(0)

names                = ['RA', 'DEC', 'ZSURV', 'ZOBS', 'MRH', 'RMAG_DRED', 'REFGMR0P1', 'GTYPE', 'HMASS', 'NMOCK', 'HALOID']
data                 = Table(np.c_[ra[isin], dec[isin], z[isin], zobs[isin], M[isin], r[isin], gmr[isin], tt[isin], hmass[isin], mxxlid[isin]], names=names) 
data['HALOID']       = data['HALOID'].data.astype(np.int) 
data['LUMDIST']      = cosmo.luminosity_distance(dat['ZSURV'].data)
data['DISTMOD']      = distmod(dat['ZSURV'].data)
data['FIELD']        = gama_field(data['RA'], data['DEC'])
data['IN_D8LUMFN']   = np.zeros_like(data['FIELD'], dtype=int)
data['CONSERVATIVE'] = np.zeros_like(data['FIELD'], dtype=int)

xyz                  = cartesian(data['RA'], data['DEC'], dat['ZSURV'])

data['CARTESIAN_X']  = xyz[:,0]
data['CARTESIAN_Y']  = xyz[:,1]
data['CARTESIAN_Z']  = xyz[:,2]

xyz                  = rotate(data['RA'], data['DEC'], xyz)

data['ROTCARTESIAN_X'] = xyz[:,0]
data['ROTCARTESIAN_Y'] = xyz[:,1]
data['ROTCARTESIAN_Z'] = xyz[:,2]

# data['GMR']          = dat['GMAG_DRED_SDSS'] - dat['RMAG_DRED_SDSS']
data['DETMAG']         = dat['RMAG_DRED']

isin, idx        = is_point_in_desi(tiles, data['RA'].data, data['DEC'].data, return_tile_index=True)

data['TILEID']   = tiles['TILEID'].data[idx]
data['ROSETTE']  = tiles['ROSETTE'].data[idx]
data['TARGETID'] = np.arange(len(data), dtype=np.int64)

root  = "/global/cscratch1/sd/mjwilson/desi/BGS/lumfn/MXXL/"
fpath = root + "galaxy_catalogue_sv3s_v4.fits"

data.write(fpath, format='fits', overwrite=True)

print('Done {}.'.format(fpath))
