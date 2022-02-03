import os
import runtime
import numpy as np
import astropy.units as u

from   astropy.coordinates import SkyCoord
from   astropy.table import Table, vstack, hstack, unique, join
from   ros_tools import tile2rosette, calc_rosr
from   gama_limits import gama_field
from   desiutil.dust import mwdust_transmission
from   desitarget.sv3.sv3_targetmask import desi_mask, bgs_mask
from   cartesian import cartesian, rotate
from   cosmo           import cosmo, distmod


root   = os.environ['DESI_ROOT'] + '/spectro/redux/everest/healpix/'
tpix   = Table.read(root + 'tilepix.fits')

tiles  = np.arange(1000)
ros    = np.array([tile2rosette(x) for x in tiles])

# https://desi.lbl.gov/trac/wiki/SurveyOps/OnePercent
# G12: [1,2]; G15: [8,9,10, 17]

uros   = np.unique(ros)
uros   = uros[uros>-1]

gama   = np.isin(ros, uros)

tiles  = tiles[gama]

tpix   = tpix[np.isin(tpix['TILEID'].data, tiles)]
hps    = np.unique(tpix['HEALPIX'].data)

root  += '/sv3/bright/'

fpaths = [root + '{}/{}/redrock-sv3-bright-{}.fits'.format(str(x)[:3], x, x) for x in hps]
fpaths = [x for x in fpaths if os.path.exists(x)]

# e.g. 280/28027/redrock-sv3-bright-28027.fits
tabs   = []

print('Gathering DESI zs.')

for x in fpaths:
    zbest = Table.read(x, hdu='REDSHIFTS')
    fmap  = Table.read(x, hdu='FIBERMAP')
    efmap = Table.read(x, hdu='EXP_FIBERMAP')

    tids  = np.unique(zbest['TARGETID'])

    # row ordered.                                                                                                                                                                                   
    assert  np.all(zbest['TARGETID'] == fmap['TARGETID'])
    
    efmap['ROS'] = tile2rosette(efmap['TILEID'].data)
    
    efmap_tid = efmap['TARGETID', 'ROS'] 
    efmap_tid = unique(efmap_tid, keys=['TARGETID'])
    efmap_tid.sort('TARGETID')

    assert  np.all(efmap_tid['TARGETID'].data == tids)
    
    del  fmap['TARGETID']
    
    zbest = hstack([zbest, fmap])
    
    # 
    zbest = join(zbest, efmap_tid, join_type='left', keys='TARGETID')
    zbest['ROS_DIST'] = 1.e4 * np.ones_like(zbest['Z'])
    
    for rosn in np.unique(efmap['ROS'].data):
        # TODO: Small rosette overlap.                                                                                                                                                               
        isin  = zbest['ROS'].data == rosn
        zbest['ROS_DIST'][isin] = calc_rosr(rosn, zbest['TARGET_RA'].data[isin], zbest['TARGET_DEC'].data[isin])
    
    tabs.append(zbest)
    
desi_zs = vstack(tabs)

# remove skies.                                                                                                                                                                                   
print('Sky frac: 1-{:.6f}'.format(np.mean(desi_zs['TARGETID'].data >= 0)))

desi_zs = desi_zs[desi_zs['TARGETID'].data >= 0]    
desi_zs = desi_zs[(desi_zs['SV3_BGS_TARGET'].data & bgs_mask['BGS_BRIGHT']) != 0]

desi_zs['GOOD_Z']     = (desi_zs['ZWARN'] == 0) & (desi_zs['DELTACHI2'] > 40)

desi_zs['ZDESI'] = desi_zs['Z']

del  desi_zs['Z']

##  Cut DESI to good redshifts.                                                                                                                                                                    
desi_zs['GAMA_FIELD'] = gama_field(desi_zs['TARGET_RA'].data, desi_zs['TARGET_DEC'].data)

archive = Table(desi_zs, copy=True)

xyz     = cartesian(desi_zs['TARGET_RA'].data, desi_zs['TARGET_DEC'].data, desi_zs['ZDESI'].data)

desi_zs['CARTESIAN_X'] = xyz[:,0]
desi_zs['CARTESIAN_Y'] = xyz[:,1]
desi_zs['CARTESIAN_Z'] = xyz[:,2]

xyz     = rotate(desi_zs['TARGET_RA'].data, desi_zs['TARGET_DEC'].data, xyz)

desi_zs['ROTCARTESIAN_X'] = xyz[:,0]
desi_zs['ROTCARTESIAN_Y'] = xyz[:,1]
desi_zs['ROTCARTESIAN_Z'] = xyz[:,2]

desi_zs['LUMDIST'] = cosmo.luminosity_distance(desi_zs['ZDESI'].data)
desi_zs['DISTMOD'] = distmod(desi_zs['ZDESI'].data)

##  HACK: PHOTSYS ASSUMED S. 
desi_zs['RMAG_DRED']  = 22.5 - 2.5 * np.log10(desi_zs['FLUX_R'].data / mwdust_transmission(desi_zs['EBV'].data, 'R', 'S', match_legacy_surveys=True))
desi_zs['IN_GOLD']    = desi_zs['GOOD_Z'].data & (desi_zs['ZDESI'] > 0.039)  & (desi_zs['ZDESI'] < 0.263)

desi_zs.pprint()

root  = os.environ['CSCRATCH'] + '/norberg/'
fpath = root + '/GAMA4/gama_gold.fits'

opath = fpath.replace('gama_gold', 'desi_sv3_gold')

print('Writing {}'.format(opath))

desi_zs.write(opath, format='fits', overwrite=True)

##  ----  GAMA GOLD
gold  = Table.read(fpath)
gold.pprint()

# DESI
uros          = [1,2,8,9,10,17] 
desi_zs       = desi_zs[np.isin(desi_zs['ROS'], uros)]
archive       = archive[np.isin(archive['ROS'], uros)]

c             = SkyCoord(ra=desi_zs['TARGET_RA']*u.degree, dec=desi_zs['TARGET_DEC']*u.degree)

print('Matching DESI to GAMA Gold.')

# GAMA
catalog       = SkyCoord(ra=gold['RA'], dec=gold['DEC'])
idx, d2d, d3d = c.match_to_catalog_3d(catalog)

# Now idx are indices into catalog that are the closest objects to each of the coordinates in c, d2d are the on-sky distances between them
gold_match   = gold[idx]

del  gold_match['FIELD']

desi_zs      = hstack([desi_zs, gold_match])
desi_zs['GAMA_SEP'] = d2d.to(u.arcsec)
desi_zs['GOOD_MATCH'] = desi_zs['GAMA_SEP'] < 0.3

max_sep      = 0.5 * u.arcsec

print('Fraction desi matched to gold at 0.5 arcseconds: {:.6f}'.format(np.mean(desi_zs['GAMA_SEP'] < max_sep)))

opath = fpath.replace('gama_gold', 'desi_gama')

print('Writing {}'.format(opath))

desi_zs.write(opath, format='fits', overwrite=True)

## --------------------
idx, d2d, d3d = catalog.match_to_catalog_3d(c)

desi_match    = archive[idx]

gold          = hstack([gold, archive])
gold['DESI_SEP'] = d2d.to(u.arcsec)
gold['GOOD_MATCH'] = gold['DESI_SEP'] < 0.3

gold['ROS_DIST'] = 1.e99

f<or rosn in uros:
    new_dist = calc_rosr(rosn, gold['RA'].data, gold['DEC'].data)
    
    gold['ROS_DIST'] = np.minimum(gold['ROS_DIST'].data, new_dist)

max_sep      = 0.5 * u.arcsec

print('Fraction gold matched to desi at 0.5 arcseconds: {:.6f}'.format(np.mean(gold['DESI_SEP'] < max_sep)))

opath = fpath.replace('gama_gold', 'gama_desi')

print('Writing {}'.format(opath))

gold.write(opath, format='fits', overwrite=True)

## --------------------
desi_zs = desi_zs[desi_zs['IN_GOLD']]

opath  = fpath.replace('gama', 'desi')

print('Writing {}'.format(opath))

desi_zs.write(opath, format='fits', overwrite=True)

print('Done.')
