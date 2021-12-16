import numpy             as np
import matplotlib.pyplot as plt
import cosmo             as cosmo
import os

from smith_kcorr  import GAMA_KCorrection
from rest_gmr     import smith_rest_gmr
from tmr_ecorr    import tmr_ecorr, tmr_q
from abs_mag      import abs_mag

kcorr_r = GAMA_KCorrection(band='R')
kcorr_g = GAMA_KCorrection(band='G')

gmr_values = np.array([0.158, 0.298, 0.419, 0.553, 0.708, 0.796, 0.960])
lim_values = [12, 19.8]
Q_values = ['RED', 'BLUE', 'ALL']
zs = np.arange(0.01, 0.6, 0.01)
mus = cosmo.distmod(zs)
aall = True

root  = os.environ['CSCRATCH'] + '/norberg/'

for k in range(len(lim_values)):
    rs = np.ones(len(zs)) * lim_values[k]

    for j in range(len(Q_values)):
        if j == 2:
            aall = True
        else:
            aall = False
            
        for idx in range(len(gmr_values)):
            gmr_0P1  = np.ones(len(zs))* gmr_values[idx]             
            gmr_0P0 = gmr_0P1 - (kcorr_g.k_nonnative_zref(0.0, zs, gmr_0P1) - kcorr_r.k_nonnative_zref(0.0, zs, gmr_0P1))
            ks = kcorr_r.k_nonnative_zref(0.0, zs, gmr_0P1)
            es = tmr_ecorr(zs, gmr_0P0, aall=aall)
            Mrs = abs_mag(rs, mus, ks, es)
            
            if (idx == 0) and (j == 0) and (k == 0):
                plt_colour = 'red'
                alpha = 1
                lw = 1
                print(gmr_values[idx], Q_values[j], lim_values[k])
                
                dat = Table()
                dat['Z'] = zs
                dat['MR_0P0'] = Mrs
                
                gmr_0P1 = "{0:.3g}".format(gmr_0P1[0])
                gmr_0P0 = "{0:.3g}".format(gmr_0P0[0])

                dat.write(root + 'GAMA4/ddrp_limits/ddrp_limit_{}_{}.fits'.format(gmr_0P1, gmr_0P0), format='fits', overwrite=True)
                
            elif (idx == 6) and (j == 2) and (k == 1):
                plt_colour = 'blue'
                alpha = 1
                lw = 1
                print(gmr_values[idx], Q_values[j], lim_values[k])
                
                dat = Table()
                dat['Z'] = zs
                dat['MR_0P0'] = Mrs
                
                gmr_0P1 = "{0:.3g}".format(gmr_0P1[0])
                gmr_0P0 = "{0:.3g}".format(gmr_0P0[0])
            
                dat.write(root + 'GAMA4/ddrp_limits/ddrp_limit_{}_{}.fits'.format(gmr_0P1, gmr_0P0), format='fits', overwrite=True)
            
            else:
                plt_colour = 'black'
                alpha = 1
                lw = 0.25

            plt.plot(zs, Mrs, color=plt_colour, lw=lw, alpha=alpha, label='gmr={}, Q={}, r={}'.format(gmr_values[idx], Q_values[j], lim_values[k]))

            
plt.xlabel('Z')
plt.ylabel(r'$M_r$ - 5 log(h)', fontsize=16)
plt.gca().invert_yaxis()
plt.show()