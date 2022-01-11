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
else:
    lumfn = lumfn[lumfn['PHI_N'] > 0.]

lumfn.pprint()

def chi2(log10phistar, Mstar, alpha):
    phistar = 10.**log10phistar
    
    res  = (lumfn['PHI_IVMAX'] - schechter(lumfn['MEDIAN_M'], phistar, Mstar, alpha))
    res /= lumfn['PHI_IVMAX_ERROR']
    res  = res * res

    return np.sum(res)
    
def lnlike(log10phistar, Mstar, alpha):
    return -chi2(log10phistar, Mstar, alpha) / 2.

x2 = chi2(-2.01, -20.89, -1.25)
print(x2)

info = {"likelihood": {"schechter": lnlike}}

info["params"] = {
     #"log10phistar": {"prior": {"min": -2.5, "max": 0.0},   "ref": -2.0,   "proposal": 0.01},  # TMR ref.  -2.01
     #"Mstar":        {"prior": {"min": -21., "max": -20.5}, "ref": -20.75, "proposal": 0.01},  # TMR ref.  -20.89 
     #"alpha":        {"prior": {"min": -1.3, "max": -1.2},  "ref": -1.25,  "proposal": 0.01}}  # TMR ref.  -1.25

     "log10phistar": {"prior": {"dist": "norm", "loc": -2.00, "scale": 0.25},   "ref": -2.00,   "proposal": 0.01},  # TMR ref.  -2.01
     "Mstar":        {"prior": {"dist": "norm", "loc": -20.89, "scale": 0.15}, "ref": -20.89, "proposal": 0.01},  # TMR ref.  -20.89 
     "alpha":        {"prior": {"dist": "norm", "loc": -1.25, "scale": 0.05},  "ref": -1.25,  "proposal": 0.01}}  # TMR ref.  -1.25

info["sampler"] = {"mcmc": {"Rminus1_stop": 0.001, "max_tries": 5000}}

updated_info, sampler = run(info, output='{}/cobaya/schechter_chain'.format(root))

print('Written to {}/cobaya/schechter_chain*'.format(root))