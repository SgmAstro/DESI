import os
import sys
import numpy as np
import argparse

from   astropy.io    import ascii
from   astropy.table import Table
from   ddp           import tmr_DDP1, tmr_DDP2, tmr_DDP3


home  = os.environ['HOME']

sys.path.append('{}/DESI'.format(home))

from   delta8_limits import delta8_tier, d8_limits

dat      = Table.read(os.environ['GOLD_DIR'] + 'gama_gold_ddp.fits')
names    = ['ZMIN', 'ZMAX', 'DDP1ZLIMS_NGAL', 'VZ', 'DENS']

tmr_DDPs = np.array([tmr_DDP1, tmr_DDP2, tmr_DDP3])

result   = Table()
rows     = []

print('\n\n')

for ddp in np.arange(1, 4, 1):
    row = [ddp, tmr_DDPs[ddp-1][0], tmr_DDPs[ddp-1][1]]
    
    for col in names:
        row += [dat.meta['DDP{}_{}'.format(ddp, col)]]

    row = tuple(row)                
    rows.append(row)

# TODO:  Note NGAL in TMR is within DDP redshift limits, not NGAL to DDP.   
names  = ['DDP', 'MIN_M', 'MAX_M'] + names
result = Table(rows=rows, names=names)
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
