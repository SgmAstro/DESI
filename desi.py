import os
import runtime
import numpy as np
import astropy.units as u

from   astropy.coordinates import SkyCoord
from   astropy.table import Table, vstack, hstack, unique, join
from   ros_tools import tile2rosette, calc_rosr
from   gama_limits import gama_field


root   = os.environ['DESI_ROOT'] + '/spectro/redux/everest/healpix/'
tpix   = Table.read(root + 'tilepix.fits')

tiles  = np.arange(1000)
ros    = np.array([tile2rosette(x) for x in tiles])

# https://desi.lbl.gov/trac/wiki/SurveyOps/OnePercent
# G12: [1,2]; G15: [8,9,10, 17]

gama   = np.isin(ros, [1,2,8,9,10,17])

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
    
    efmap['ROS']      = tile2rosette(efmap['TILEID'].data)
    
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

desi_zs.pprint()

## GOLD:
root = os.environ['CSCRATCH'] + '/norberg/'
fpath = root + '/GAMA4/gama_gold.fits'

gold = Table.read(fpath)
gold.pprint()

##  Cut DESI to good redshifts.
desi_zs       = desi_zs[(desi_zs['ZWARN'] == 0) & (desi_zs['DELTACHI2'] > 40)]
desi_zs       = desi_zs[(desi_zs['Z'] > 0.039)  & (desi_zs['Z'] < 0.263)]

desi_zs['GAMA_FIELD'] = gama_field(desi_zs['TARGET_RA'], desi_zs['TARGET_DEC'])

# DESI
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

max_sep      = 0.5 * u.arcsec

print('Fraction matched to 0.5 arcseconds: {:.6f}'.format(np.mean(desi_zs['GAMA_SEP'] < max_sep)))

opath = fpath.replace('gama', 'desi')

print('Writing {}'.format(opath))

desi_zs.write(opath, format='fits', overwrite=True)

print('Done.')
