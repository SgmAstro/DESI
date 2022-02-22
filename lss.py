import fitsio
import numpy as np

from   astropy.table import Table


def fetch_lss(version=2.1, pprint=False, sort=False):
    # https://desi.lbl.gov/trac/wiki/ClusteringWG/LSScat/SV3/version2.1/
    #
    # /global/cfs/cdirs/desi/survey/catalogs/SV3/LSS/everest/LSScats/2.1/

    # ('RA', 'DEC', 'TARGETID', 'Z', 'NTILE', 'TILES', 'rosette_number', 'rosette_r', 'FRACZ_TILELOCID', 'BITWEIGHTS', 'PROB_OBS', 'WEIGHT_ZFAIL', 'WEIGHT', 'flux_r_dered', 'NZ')
    clustering = Table.read(f'/global/cfs/cdirs/desi/survey/catalogs/SV3/LSS/everest/LSScats/{version}/BGS_BRIGHT_clustering.dat.fits')

    # NTILE, TILES, TILELOCIDS, LOCATION_ASSIGNED, TILELOCID_ASSIGNED, sort, COMP_TILE, rosette_number, rosette_r, FRACZ_TILELOCID, BITWEIGHTS, PROB_OBS
    full       = Table.read(f'/global/cfs/cdirs/desi/survey/catalogs/SV3/LSS/everest/LSScats/{version}/BGS_BRIGHT_full.dat.fits')

    if pprint:
        print(clustering.dtype.names)
        
        clustering.pprint()
        kp3_cat.pprint()

        full_cols = np.array(full.dtype.names).astype(str)

        print(full_cols)

        # lumfn.io
        # np.savetxt('full_cols.txt', full_cols, fmt='%s')

    if sort:
        clustering.sort('TARGETID')
        full.sort('TARGETID')
        
    return  clustering, full


if __name__ == '__main__':
    fetch_lss()
    
