import numpy             as np
import pylab             as pl
import matplotlib.pyplot as plt

from   speclite.filters import load_filters
from   scipy.interpolate import interp1d


def central(wave, response):
    return np.sum(wave * response) / np.sum(response)

zref    = 0.1

bass_g  = np.loadtxt('bass-g.txt')
bass_r  = np.loadtxt('bass-r.txt')

decam_g = np.loadtxt('decam.g.am1p4.dat.txt')
decam_r = np.loadtxt('decam.r.am1p4.dat.txt')
decam_z = np.loadtxt('decam.z.am1p4.dat.txt')

sdss    = load_filters('sdss2010-*')
sdss_g  = np.c_[sdss[1].wavelength, sdss[1].response]
sdss_r  = np.c_[sdss[2].wavelength, sdss[2].response]
sdss_z  = np.c_[sdss[4].wavelength, sdss[4].response]

fs      = [bass_g, bass_r, decam_g, decam_r, decam_z, sdss_g, sdss_r, sdss_z]
ls      = ['bass g', 'bass r', 'decam g', 'decam r', 'decam z', 'sdss g', 'sdss r', 'sdss z']    
cs      = plt.rcParams['axes.prop_cycle'].by_key()['color']
    
pl.plot(bass_g[:,0], bass_g[:,-1], label=ls[0],   c=cs[0])
pl.plot(bass_r[:,0], bass_r[:,-1], label=ls[1],   c=cs[1])

pl.plot(decam_g[:,0], decam_g[:,-1], label=ls[2], c=cs[2])
pl.plot(decam_r[:,0], decam_r[:,-1], label=ls[3], c=cs[3])
pl.plot(decam_z[:,0], decam_z[:,-1], label=ls[4], c=cs[4])

pl.plot(sdss_g[:,0], sdss_g[:,-1], label=ls[5], c=cs[5])
pl.plot(sdss_r[:,0], sdss_r[:,-1], label=ls[6], c=cs[6])
pl.plot(sdss_z[:,0], sdss_z[:,-1], label=ls[7], c=cs[7])

baseline = central(decam_r[:,0], decam_r[:,-1])

stack = {}

for ll, ff, cc in zip(ls, fs, cs):
    wave = ff[:,0]
    response = ff[:,-1]

    median = central(wave, response)

    # ref. frame freq. = beta * obs freq. 
    
    beta = baseline / median
    
    pl.axvline(median, lw=1., color=cc, label=r'$\beta$: {:.2f}'.format(beta))    

    # print(ll.ljust(15), zref, median, '{:.4f}'.format(beta * (1. + zref) - 1.))

    stack[ll] = {'central': median, 'wave': wave, 'response': response}
    
# pl.legend(frameon=False, ncol=2)
# pl.show()

# print(stack.keys())

def similarity(obs_filter='decam r', rest_filter='sdss r', ttype='correlation'):
    ## --  Cross-correlate -- 
    obs = stack[obs_filter]
    rest = stack[rest_filter]

    obs = interp1d(obs['wave'] - obs['central'], obs['response'], kind='linear', axis=-1, copy=True, bounds_error=False, fill_value=0.0, assume_sorted=False)
    rest = interp1d(rest['wave'] - rest['central'], rest['response'], kind='linear', axis=-1, copy=True, bounds_error=False, fill_value=0.0, assume_sorted=False)

    wave = np.arange(-1.e4, 1.e4, 1.)

    obs = obs(wave)
    rest = rest(wave)
    
    # https://numpy.org/doc/stable/reference/generated/numpy.correlate.html
    result  = np.sum(obs * rest) / len(obs)
    result /= np.sqrt(np.sum(obs**2.) / len(obs))
    result /= np.sqrt(np.sum(rest**2.) / len(rest)) 

    return result

f = open('filters.txt', 'w')

f.write('# obs_filter\trest_filter\tobs_central\trest_central\tbeta\t\tsimilarity\n')

for l1 in ls:
    for l2 in ls:
        rr = similarity(l1, l2)

        beta = stack[l1]['central'] / stack[l2]['central']
        
        f.write('{}\t{}\t{:.6e}\t{:.6e}\t{:.6e}\t{:.6e}\n'.format(l1.ljust(15), l2.ljust(15), stack[l1]['central'], stack[l2]['central'], beta, rr))
        
f.close()
