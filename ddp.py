import os
import fitsio
import numpy             as np
import matplotlib.pyplot as plt

from   cosmo             import volcom
from   scipy.interpolate import interp1d


tmr_DDP1       = [-21.8, -20.1]
tmr_DDP2       = [-20.6, -19.3]
tmr_DDP3       = [-19.6, -17.8]

root           = os.environ['GOLD_DIR'] + '/ddrp_limits/'

_bright_curve  = fitsio.read(root + '/ddrp_limit_3.fits')  #  7 (12.0 QCOLOR 0.131)
_faint_curve   = fitsio.read(root + '/ddrp_limit_17.fits') # 27 (19.8 QCOLOR 1.067)

# TODO: extend the curve limits and put bounds_error back on.
bright_curve   = interp1d(_bright_curve['M0P0_QALL'], _bright_curve['Z'], kind='linear', copy=True, bounds_error=False, fill_value=0.0, assume_sorted=False)
bright_curve_r = interp1d(_bright_curve['Z'],         _bright_curve['M0P0_QALL'], kind='linear', copy=True, bounds_error=False, fill_value=0.0, assume_sorted=False)

faint_curve    = interp1d(_faint_curve['M0P0_QALL'],  _faint_curve['Z'],  kind='linear', copy=True, bounds_error=False, fill_value=1.0, assume_sorted=False)
faint_curve_r  = interp1d(_faint_curve['Z'],          _faint_curve['M0P0_QALL'],   kind='linear', copy=True, bounds_error=False, fill_value=1.0, assume_sorted=False)

def get_ddps(Area, M_0P0s, zs):
    result   = np.zeros(len(zs) * 3, dtype=int).reshape(len(zs), 3)
    zlims    = {}
    
    for i, lims in enumerate([tmr_DDP1, tmr_DDP2, tmr_DDP3]):
        in_ddp  = (M_0P0s >= lims[0]) & (M_0P0s <= lims[1])

        zmax    = np.atleast_1d(faint_curve(lims[1]))[0]
        zmin    = np.atleast_1d(bright_curve(lims[0]))[0]
        
        exclude = (zs > zmax) | (zs < zmin)

        in_ddp  = in_ddp & ~exclude  
        in_ddpz = ~exclude
        
        result[in_ddp, i] = 1

        ddp_zs  = zs[in_ddp]

        # print(zmin, zmax, len(ddp_zs))
        
        zmax = np.array([zmax, ddp_zs.max()]).min()
        zmin = np.array([zmin, ddp_zs.min()]).max()
        
        zlims['DDP{}_ZMIN'.format(i+1)] = zmin
        zlims['DDP{}_ZMAX'.format(i+1)] = zmax

        zlims['DDP{}_VZ'.format(i+1)]   = volcom(zmax, Area) - volcom(zmin, Area)

        zlims['DDP{}ZLIMS_NGAL'.format(i+1)] = np.count_nonzero(in_ddpz)

        zlims['DDP{}_NGAL'.format(i+1)] = np.count_nonzero(in_ddp) 
        zlims['DDP{}_DENS'.format(i+1)] = np.count_nonzero(in_ddp) / zlims['DDP{}_VZ'.format(i+1)] 
                
    return  result, zlims, faint_curve_r(zs)


if __name__ == '__main__':
    print('Done.')
