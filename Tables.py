
import os
import argparse
from astropy.table import Table

parser = argparse.ArgumentParser(description='Generate Gold luminosity function.')
parser.add_argument('-f', '--field', type=str, help='Select equatorial GAMA field: G9, G12, G15', required=False, default=None)
    
user = os.environ['USER']

# TODO: Add to headers.
tmr_DDP1     = [-21.8, -20.1]
tmr_DDP2     = [-20.6, -19.3]
tmr_DDP3     = [-19.6, -17.8]

Mr = np.array([tmr_DDP1, tmr_DDP2, tmr_DDP3])

rpath = '/cosma/home/durham/{}/data/GAMA4/randoms/randoms_bd_ddp_n8_{}_0.fits'.format(user, field)
rand = Table.read(rpath)

abs_mag = []
ddp_min = []
ddp_max = []
ddp_ngal = []
ddp_vol = []
ddp_dens = []
ddp = [1,2,3]

for idx in ddp:
    abs_mag.append(Mr[idx-1])
    ddp_min.append(rand.meta['DDP{}_ZMIN'.format(idx)])
    ddp_max.append(rand.meta['DDP{}_ZMAX'.format(idx)])
    ddp_ngal.append(rand.meta['DDP{}_NGAL'.format(idx)])
    ddp_vol.append(rand.meta['DDP{}_VZ'.format(idx)])
    ddp_dens.append(rand.meta['DDP{}_DENS'.format(idx)])

# Generate Table 2 of McNaught-Roberts (2014)
t1 = Table([ddp, abs_mag, ddp_min, ddp_max, ddp_ngal, ddp_vol, ddp_dens], names=('DDP', 'MRH', 'ZMIN', 'ZMAX', 'N_DDP','VOL_DDP', 'DENS_DDP'))

fscale = []
n_ddp = []
label = ['d0', 'd1', 'd2', 'd3']

print('bin', 'overdensities', 'f_delta', 'n_gal/10^3')
for idx in range(4):
    
    fpath = '/cosma/home/durham/{}/data/GAMA4/gama_gold_{}_ddp_n8_d0_{}.fits'.format(user, field, idx)
    dat = Table.read(fpath)
    n_ddp.append(len(dat)/10**3)
    
    fscale.append(rand.meta['DDP1_d{}_VOLFRAC'.format(idx)])

# Generate Table 3 of McNaught-Roberts (2014)
t2 = Table([label, d8_limits, fscale, n_ddp], names=('label', 'd8_limits', 'f_delta', 'N_DDP/10^3'))
    
print(t1)
print(t2)