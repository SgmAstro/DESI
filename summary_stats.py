import os
import sys
import argparse
import subprocess
import numpy as np

from   astropy.io    import ascii
from   astropy.table import Table


home  = os.environ['HOME']

sys.path.append('{}/DESI'.format(home))

from   ddp           import tmr_DDP1, tmr_DDP2, tmr_DDP3
from   delta8_limits import delta8_tier, d8_limits


dat              = Table.read(os.environ['GOLD_DIR'] + '/gama_gold_ddp.fits')
names            = ['ZMIN', 'ZMAX', 'VZ', 'DENS']

tmr_DDPs         = np.array([tmr_DDP1, tmr_DDP2, tmr_DDP3])
tmr_ngold_ratios = np.array([1.002, 0.591, 0.097])
tmr_ddp_ratios   = np.array([1., 0.589, 0.097])

result   = Table()
rows     = []

print('\n\n')

for ddp, tmr_DDP, tmr_ngold_ratio, tmr_ddp_ratio in zip(np.arange(1, 4, 1), tmr_DDPs, tmr_ngold_ratios, tmr_ddp_ratios):
    row = [ddp, tmr_DDP[0], tmr_DDP[1]]
    
    for col in names:
        row += [dat.meta['DDP{}_{}'.format(ddp, col)]]
        
    row += [dat.meta['DDP{}{}'.format(ddp, 'ZLIMS_NGAL')]]
    row += [1. * dat.meta['DDP{}_NGAL'.format(ddp)] / dat.meta['GOLD_NGAL']]
    row += [tmr_ngold_ratio]
        
    row = tuple(row)         
    rows.append(row)

names  = ['DDP', 'MIN_M', 'MAX_M'] + names + ['ZLIMS_NGAL', 'N_NGOLD', 'N_NREF_TMR']
result = Table(rows=rows, names=names)

result['ZLIMS_NGAL']      = result['ZLIMS_NGAL'] / 10**3
result['VZ']              = result['VZ'] / 10**6
result['DENS']            = result['DENS'] / 10**-3
result['N_GAL/N_GAL_MAX'] = 1. * result['ZLIMS_NGAL'] / max(result['ZLIMS_NGAL'])

for name in result.dtype.names:
    result[name]          = np.round(result[name], 3)
    
result.rename_column('MIN_M', '$M_{Min.}$')
result.rename_column('MAX_M', '$M_{Max.}$')
result.rename_column('ZLIMS_NGAL', '$N_{GAL} / 10^3$')
result.rename_column('VZ', '$V_{DDP}$ / 10^6$')
result.rename_column('DENS', '$\\rho_{DDP} / 10^{-3}$')
result.rename_column('N_NGOLD', '$N/N_{GOLD}$')
result.rename_column('N_NREF_TMR', '$N/N_{REF}_{TMR}$')
result.rename_column('N_GAL/N_GAL_MAX', '$N_{GAL}/N_{GALMAX}$')

result.pprint()

# https://arxiv.org/pdf/1409.4681.pdf
ascii.write(result, 'tables/Tab2.tex', Writer=ascii.Latex, latexdict=ascii.latex.latexdicts['AA'], overwrite=True)

exit(0)

rows = []

rand = Table.read('{}/randoms_bd_ddp_n8_G9_0.fits'.format(os.environ['RANDOMS_DIR']))

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

# subprocess.run('cd {}/tables; pdflatex compile_tables.tex'.format(os.environ['CODE_ROOT']), shell=True, check=True)

print('\n\nDone.\n\n')
