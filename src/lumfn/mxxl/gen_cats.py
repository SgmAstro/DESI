import os
import sys
import time
import h5py
import fitsio
import pylab             as     pl
import numpy             as     np
import astropy.io.fits   as     fits

from   astropy.table     import Table, join, vstack
from   desitarget.sv3    import sv3_targetmask
from   matplotlib.pyplot import figure
from   scipy.optimize    import curve_fit

sys.path.append('/global/homes/m/mjwilson/desi/BGS/lumfn/bin')
sys.path.append('/global/homes/m/mjwilson/desi/BGS/lumfn/py/')
sys.path.append('/global/homes/m/mjwilson/desi/BGS/lumfn/py/lumfn/')

sys.path.append('/global/homes/m/mjwilson/desi/LSS/bin/')
sys.path.append('/global/homes/m/mjwilson/desi/LSS/py/')

from   distances          import comoving_distance, comoving_volume, dist_mod
from   ajs_kcorr          import ajs_kcorr
from   mxxl_ecorr         import mxxl_ecorr
from   tmr_ecorr          import tmr_ecorr
from   abs_mag            import abs_mag
from   vmax               import zmax, vmax
from   ref_gmr            import ajs_reference_gmr, tmr_reference_gmr
from   params             import params
from   zmin               import zmin
from   rlim               import rlim
from   MXXL.define_sample import define_sample
from   MXXL.mxxl_params   import mxxl_params


root  = "/global/cscratch1/sd/mjwilson/desi/BGS/lumfn/MXXL/"

# fpath = root + "galaxy_catalogue_small.hdf5"
fpath = root + "galaxy_catalogue_sv3s_v4.fits"

print(fpath)

version          = 0.8
todisk           = True
dryrun           = False
tmr              = False
runtime_lim      = 10. # only if dryrun.

odir             = os.environ['CSCRATCH'] + '/desi/BGS/lumfn/'

print(version, odir)
'''
f   = h5py.File(fpath, mode='r')

ra  = f["Data/ra"][...]
dec = f["Data/dec"][...]
z   = f["Data/z_cos"][...]
r   = f["Data/app_mag"][...]

# Absolute magnitude: k-corrected to the r at z=0.1 band but no E.
# Rest frame g-r colour in the z=0.1 shifted SDSS bands.
Mrh = f["Data/abs_mag"][...]
gmr = f["Data/g_r"][...]

f.close()


mxxl             = Table(np.c_[ra, dec, z, r, Mrh, gmr], names=['RA', 'DEC', 'Z', 'RMAG_DRED', 'MRH', 'REFGMR0P1'])
'''

mxxl             = Table.read(fpath)
mxxl             = mxxl[np.isin(mxxl['NMOCK'].data, np.arange(0, 5, 1))]

mxxl             = mxxl[define_sample(mxxl)]

if todisk:
    mxxl.write('{}/MXXL/bright_v{:.1f}.fits'.format(odir, version), format='fits', overwrite=True)

exit(0)
    
kcorrector   = ajs_kcorr()

derived      = []
start        = time.time()

tids         = []

for ii, row in enumerate(mxxl):
    tids.append(row['TARGETID'])
    
    rmag     = row['RMAG_DRED']
    zz       = row['Z']
    
    # Needed to define r mag. limit (at 19.5 for S).
    psys     = 'S'

    insample = define_sample(row)
    
    #  See dedicated fsky calc. notebook.
    vol      = comoving_volume(zmin(params['vmin']), zz, fsky=mxxl_params['fsky'])

    # one_reference_gmr(kcorrector, gmr, zz, zref=kcorrector.z0, ecorr=False)                                                                                                                                                                
    org_gmr = row['REFGMR0P1']
    
    gmr     = org_gmr + kcorrector.ref_eval(org_gmr, zz, band='g')[0] - kcorrector.ref_eval(org_gmr, zz, band='r')[0]
    
    Mrh     = abs_mag(kcorrector, rmag, None, zz, ref_gmr=org_gmr, tmr=tmr).item()

    maxz    = zmax(kcorrector, rlim(psys), Mrh, None, zz, ref_gmr=org_gmr, tmr=tmr)

    ##  (!!!)
    maxz    = np.minimum(maxz, params['zmax'])
    
    maxv    = vmax(kcorrector, rlim(psys), Mrh, None, zz, min_z=zmin(params['vmin']), fsky=mxxl_params['fsky'], max_z=maxz)

    mu      = dist_mod(maxz)
    rk      = kcorrector.ref_eval(org_gmr, maxz, band='r')[0]

    if tmr:
        rk  = kcorrector.shift_refz(rk, org_gmr, ref_z=0.0, band='r')

        tmr_ref_gmr = tmr_reference_gmr(kcorrector, org_gmr)

        rk += tmr_ecorr(maxz, tmr_ref_gmr)

    else:
        rk += mxxl_ecorr(maxz)

    rk      = rk.item()
        
    derived.append([insample, mu, vol, Mrh, rk, gmr, maxz, 1. / maxv])

    if not (ii % 100):
        runtime = (time.time()-start) / 60.

        percentage_complete = 100. * ii / len(mxxl)

        print('{:.2f} complete after {:.2f} minutes.'.format(percentage_complete, runtime))
        
        if dryrun & (runtime > runtime_lim):
            break

derived = Table(np.array(derived), names=['INSAMPLE', 'DISTMOD_ZMAX', 'VOLUME', 'MRH', 'RKCORR_ZMAX', 'GMR_DRED', 'ZMAX', 'IVMAX'])
derived['TARGETID'] = np.array(tids, dtype=np.int64)
derived.pprint(max_width=-1)

if todisk:
    derived.write('{}/MXXL/bright_derived_v{:.1f}.fits'.format(odir, version), format='fits', overwrite=True)
    
print('Done.')
