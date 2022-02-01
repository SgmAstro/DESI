import os
import sys
import numpy as np
import argparse

from   astropy.io    import ascii
from   astropy.table import Table


home  = os.environ['HOME']

sys.path.append('{}/DESI'.format(home))

from   ddp           import tmr_DDP1, tmr_DDP2, tmr_DDP3
from   delta8_limits import delta8_tier, d8_limits

dat      = Table.read(os.environ['GOLD_DIR'] + 'gama_gold_ddp.fits')
names    = ['ZMIN', 'ZMAX', 'DDP1ZLIMS_NGAL', 'VZ', 'DENS']

tmr_DDPs = np.array([tmr_DDP1, tmr_DDP2, tmr_DDP3])
tmr_ngold_ratio = np.array([1.002, 0.591, 0.097])
tmr_ddp_ratio = np.array([1., 0.589, 0.097])


result   = Table()
rows     = []

print('\n\n')

for ddp in np.arange(1, 4, 1):
    row = [ddp, tmr_DDPs[ddp-1][0], tmr_DDPs[ddp-1][1]]
    
    for col in names:
        try:
            row += [dat.meta['DDP{}_{}'.format(ddp, col)]]
        except:
            row += [dat.meta['DDP{}{}'.format(ddp, col)]]

    
    row += [dat.meta['DDP{}_NGAL'.format(ddp)]/ dat.meta['GOLD_NGAL']]
    row += [tmr_ngold_ratio[ddp-1]]
        
    row = tuple(row)                
    rows.append(row)

    dat.meta['DDP{}_NGAL'.format(ddp)]/ dat.meta['GOLD_NGAL']

names += ['N_NGOLD', 'N_NGOLD_TMR']
    
names  = ['DDP', 'MIN_M', 'MAX_M'] + names
result = Table(rows=rows, names=names)

result['ZLIMS_NGAL'] = result['ZLIMS_NGAL'] / 10**3
result['VZ'] = result['VZ'] / 10**6
result['DENS'] = result['DENS'] / 10**-3
result['N_GAL/N_GAL_MAX'] = result['ZLIMS_NGAL'] / max(result['ZLIMS_NGAL'])

for name in names:
    result[name] = np.round(result[name], 3)
    
result.rename_column('ZLIMS_NGAL', '$N_{GAL} / 10^3$')
result.rename_column('VZ', '$V_{DDP}$ / 10^6$')
result.rename_column('DENS', '$\\rho_{DDP} / 10^{-3}$')
result.rename_column('N_NGOLD', '$N/N_{GOLD}$')
result.rename_column('N_NGOLD_TMR', '$N/N_{GOLD}_{TMR}$')

result.pprint()

# https://arxiv.org/pdf/1409.4681.pdf
ascii.write(result, 'tables/Tab2.tex', Writer=ascii.Latex, latexdict=ascii.latex.latexdicts['AA'], overwrite=True)

rows = []

rand = Table.read('{}/data/GAMA4/randoms/randoms_bd_ddp_n8_G9_0.fits'.format(os.environ['HOME']))

print('\n\n')

for idx in np.arange(9):    
    nd8 = 0 
    
    for field in ['G9', 'G12', 'G15']:
        dat  = Table.read('{}/data/GAMA4/gama_gold_{}_ddp_n8_d0_{}.fits'.format(os.environ['HOME'], field, idx))
        nd8 += len(dat) / 1.e3
    
    rows.append(('d{}'.format(idx), d8_limits[idx][0], d8_limits[idx][1], nd8, rand.meta['DDP1_d{}_VOLFRAC'.format(idx)]))

print('\n\n')

# Generate Table 3 of McNaught-Roberts (2014).
result = Table(rows=rows, names=['Label', 'Min_{d8}', 'Max_{d8}', 'N_{d8} [1e3]', 'fd8'])
result.pprint()

# https://arxiv.org/pdf/1409.4681.pdf
ascii.write(result, 'tables/Tab3.tex', Writer=ascii.Latex, latexdict=ascii.latex.latexdicts['AA'], overwrite=True)

print('\n\nDone.\n\n')
