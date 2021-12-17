import  os
import  numpy as np
import  getdist.plots as gdplt

from    schechter import schechter
from    astropy.table import Table
from    cobaya.run import run
from    scipy import stats
from    getdist.mcsamples import MCSamplesFromCobaya

# Run me on interactive:
# desienv master
# srun -N 1 -n 1 python fit_schechter.py

known = False

root  = os.environ['CSCRATCH'] + '/norberg/GAMA4/'
fpath = root + '/gama_gold_lumfn.fits'
lumfn = Table.read(fpath)

if known:
    lumfn['PHI_N'] = schechter(lumfn['MEDIAN_M'], -2.01, -20.89, -1.25)
    lumfn['PHI_N_ERROR'] = 1.e-2 * lumfn['PHI_N']

lumfn.pprint()

def chi2(log10phistar, Mstar, alpha):
    phistar = 10.**log10phistar
    
    res  = (lumfn['PHI_IVMAX'] - schechter(lumfn['MEDIAN_M'], phistar, Mstar, alpha))
    res /= lumfn['PHI_IVMAX_ERROR']
    res  = res * res

    return np.sum(res)

def lnlike(log10phistar, Mstar, alpha):
    return -chi2(log10phistar, Mstar, alpha) / 2.
    
info = {"likelihood": {"schechter": lnlike}}

info["params"] = {
     "log10phistar": {"prior": {"min": -2.5, "max": 0.0},   "ref": -2.0,   "proposal": 0.01},  # TMR ref.  -2.01
     "Mstar":        {"prior": {"min": -21., "max": -20.5}, "ref": -20.75, "proposal": 0.01},  # TMR ref.  -20.89 
     "alpha":        {"prior": {"min": -1.3, "max": -1.2},  "ref": -1.25,  "proposal": 0.01}}  # TMR ref.  -1.25

info["sampler"] = {"mcmc": {"Rminus1_stop": 0.001, "max_tries": 1000}}

updated_info, sampler = run(info, output='{}/cobaya/schechter_chain'.format(root))

print('Written to {}/cobaya/schechter_chain*'.format(root))

# gdsamples = MCSamplesFromCobaya(updated_info, sampler.products()["sample"])
# gdplot    = gdplt.get_subplot_plotter(width_inch=5)

# gdplot.triangle_plot(gdsamples, ["log10phistar", "Mstar", "alpha"], filled=True)

# gdplot    = gdplt.get_subplot_plotter(width_inch=5)
# gdplot.plots_1d(gdsamples, ["r", "theta"], nx=2)
