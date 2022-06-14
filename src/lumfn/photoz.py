import pylab               as     pl
import numpy               as     np 
import astropy.units       as     u

from   cosmo               import distcom
from   astropy.table       import Table, join, unique
from   astropy.coordinates import SkyCoord


sci        = Table.read('/cosma/home/durham/dc-wils7/data/gkvScienceCatv02.fits')
sci        = sci[::10]

'''
for x in sci.dtype.names:
    print(x)
'''
# photoz   = Table.read('/cosma/home/durham/dc-wils7/data/EAZYPhotoZv02.fits')

cs             = SkyCoord(ra=sci['RAcen'], dec=sci['Deccen'])
catalogs       = SkyCoord(ra=sci['RAcen'], dec=sci['Deccen'])

for cc in cs:
    idxc, idxcatalog, d2d, d3d = catalog.search_around_sky(cs, 6. * u.arcmin)

    matches        = sci[idxc[d2d > 2. * u.arcsec]]
    matches['D2D'] = d2d[d2d > 2. * u.arcsec] 
    matches.sort('D2D')
    
    matches        = unique(sci, keys='uberID', keep='first')
    sci.pprint()

# sci    = join(sci, photoz, join_type='left', keys='uberID')
# sci    = sci[sci['SC'] >= 7]
# sci    = sci[sci['NQ'] >= 3]
# sci.pprint()

'''
chi_s  = distcom(sci['Z'])
chi_p  = distcom(sci['z_peak'])

pl.plot(sci['Z'], np.abs(chi_s - chi_p), marker=',', lw=0.0, c='k')

pl.xlim(0.0, 0.8)
# pl.ylim(0.0, 0.8)

pl.xlabel(r'$z_s$')
# pl.ylabel(r'$z_p$')

pl.yscale('log')

pl.savefig('test.pdf')
'''
