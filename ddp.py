import os
import fitsio
import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate

def ddp(M_0P0, z):
    tmr_DDP1Mr = [-20.1, -21.8]
    tmr_DDP2Mr = [-19.3, -20.6]
    tmr_DDP3Mr = [-17.8, -19.6]

    gmr_0P1, gmr_0P0 = 0.158, 0.16
    fpath = os.environ['CSCRATCH'] + '/norberg/GAMA4/ddrp_limits/ddrp_limit_{}_{}.fits'.format(gmr_0P1, gmr_0P0)
    curve2 = fitsio.read(fpath)

    gmr_0P1, gmr_0P0 = 0.96, 0.96
    fpath = os.environ['CSCRATCH'] + '/norberg/GAMA4/ddrp_limits/ddrp_limit_{}_{}.fits'.format(gmr_0P1, gmr_0P0)
    curve1 = fitsio.read(fpath)

    f1 = interpolate.interp1d(curve1['MR_0P0'], curve1['Z'], kind='linear', copy=True, bounds_error=True, fill_value=np.NaN, assume_sorted=False)

    f2 = interpolate.interp1d(curve2['MR_0P0'], curve2['Z'], kind='linear', copy=True, bounds_error=True, fill_value=np.NaN, assume_sorted=False)

    Mr_new1 = np.arange(curve1['MR_0P0'].min(), curve1['MR_0P0'].max(), 0.5)
    Mr_new2 = np.arange(curve2['MR_0P0'].min(), curve2['MR_0P0'].max(), 0.5)

    tmr_DDP1z = [f1(tmr_DDP1Mr), f2(tmr_DDP1Mr)]
    tmr_DDP2z = [f1(tmr_DDP2Mr), f2(tmr_DDP2Mr)]
    tmr_DDP3z = [f1(tmr_DDP3Mr), f2(tmr_DDP3Mr)]

    DDP1_true = (M_0P0 < tmr_DDP1Mr[0]) & (M_0P0 > tmr_DDP1Mr[1]) & (z < tmr_DDP1z[0][0]) & (z > tmr_DDP1z[1][1])
    DDP2_true = (M_0P0 < tmr_DDP2Mr[0]) & (M_0P0 > tmr_DDP2Mr[1]) & (z < tmr_DDP2z[0][0]) & (z > tmr_DDP2z[1][1])
    DDP3_true = (M_0P0 < tmr_DDP3Mr[0]) & (M_0P0 > tmr_DDP3Mr[1]) & (z < tmr_DDP3z[0][0]) & (z > tmr_DDP3z[1][1])

    # returns [0, 1, 0] array
    # need to covert bool to int
    return [DDP1_true, DDP2_true, DDP3_true]


