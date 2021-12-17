import os
import numpy as np
import astropy.units as u

from astropy.coordinates import SkyCoord
from astropy.table import Table, vstack, hstack


root = os.environ['DESI_ROOT'] + '/spectro/redux/everest/healpix/'
tpix = Table.read(root + 'tilepix.fits')

# https://github.com/desihub/LSS/blob/3d65b099a763a41179cb8e706f85482b3633a254/py/LSS/SV3/cattools.py#L21
def tile2rosette(tile):
    if tile < 433:
        return (tile-1)//27
    else:
        if tile >= 433 and tile < 436:
            return 13
        if tile >= 436 and tile < 439:
            return 14
        if tile >= 439 and tile < 442:
            return 15
        if tile >= 442 and tile <=480:
            return (tile-442)//3
            
        if tile > 480:
            return tile//30    
    return 999999 #shouldn't be any more?


tiles = np.arange(1000)
ros   = np.array([tile2rosette(x) for x in tiles])

# https://desi.lbl.gov/trac/wiki/SurveyOps/OnePercent
# G12: [1,2]; G15: [8,9,10, 17]

gama  = np.isin(ros, [1,2,8,9,10,17])

tiles = tiles[gama]

tpix  = tpix[np.isin(tpix['TILEID'].data, tiles)]
hps   = np.unique(tpix['HEALPIX'].data)

root += '/sv3/bright/'

fpaths = [root + '{}/{}/redrock-sv3-bright-{}.fits'.format(str(x)[:3], x, x) for x in hps]
fpaths = [x for x in fpaths if os.path.exists(x)]

# e.g. 280/28027/redrock-sv3-bright-28027.fits
tabs = []

print('Gathering DESI zs.')

for x in fpaths:
    zbest = Table.read(x, hdu='REDSHIFTS')
    fmap  = Table.read(x, hdu='FIBERMAP')

    # row ordered.
    assert np.all(zbest['TARGETID'] == fmap['TARGETID'])

    del fmap['TARGETID']
    
    zbest = hstack([zbest, fmap])

    # remove skies.
    zbest = zbest[zbest['TARGETID'] >= 0]
    
    tabs.append(zbest)
    
desi_zs = vstack(tabs)
desi_zs.pprint()

## GOLD:
root = os.environ['CSCRATCH'] + '/norberg/'
fpath = root + '/GAMA4/gama_gold.fits'

gold = Table.read(fpath)
gold.pprint()

# DESI
c             = SkyCoord(ra=desi_zs['TARGET_RA']*u.degree, dec=desi_zs['TARGET_DEC']*u.degree)

print('Matching DESI to GAMA Gold.')

# GAMA
catalog       = SkyCoord(ra=gold['RA'], dec=gold['DEC'])
idx, d2d, d3d = c.match_to_catalog_3d(catalog)

# Now idx are indices into catalog that are the closest objects to each of the coordinates in c, d2d are the on-sky distances between them
gold_match   = gold[idx]

desi_zs      = hstack([desi_zs, gold_match])
desi_zs['GAMA_SEP'] = d2d

max_sep      = 1.0 * u.arcsec

print(np.mean(desi_zs['GAMA_SEP'] < max_sep))

opath = fpath.replace('gama', 'desi')

print('Writing {}'.format(opath))

desi_zs.write(opath, format='fits', overwrite=True)

print('Done.')
