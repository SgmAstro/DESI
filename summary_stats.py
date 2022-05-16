import os
import sys
import argparse
import subprocess
import numpy as np

from   astropy.io    import ascii
from   astropy.table import Table
from   ddp           import tmr_DDP1, tmr_DDP2, tmr_DDP3
from   delta8_limits import delta8_tier, d8_limits
from   findfile      import findfile

parser = argparse.ArgumentParser(description='Generate Summary Stats')
parser.add_argument('-s', '--survey', help='Select survey', default='gama')

args             = parser.parse_args()
survey           = args.survey

fpath            = findfile(ftype='ddp', survey=survey)
dat              = Table.read(fpath)
names            = ['ZMIN', 'ZMAX', 'VZ', 'DENS']

tmr_DDPs = [tmr_DDP1, tmr_DDP2, tmr_DDP3]

result = Table()
rows   = []

print('\n\n')

for ddp, tmr_DDP in zip(np.arange(1, 4, 1), tmr_DDPs):
    row = [ddp, tmr_DDP[0], tmr_DDP[1]]
    
    for col in names:
        row += [dat.meta['DDP{}_{}'.format(ddp, col)]]
        
    row += [dat.meta['DDP{}{}'.format(ddp, 'ZLIMS_NGAL')]]        
    row  = tuple(row)         
    rows.append(row)

names  = ['DDP', 'MIN_M', 'MAX_M'] + names + ['ZLIMS_NGAL']
result = Table(rows=rows, names=names)

result['ZLIMS_NGAL']      = result['ZLIMS_NGAL'] / 10**3
result['VZ']              = result['VZ'] / 10**6
result['DENS']            = result['DENS'] / 10**-3

for name in result.dtype.names:
    result[name]          = np.round(result[name], 4)

opath = 'tables/Tab2.fits'
result.write(opath, format='fits', overwrite=True)
    
result.rename_column('MIN_M', r'$M_{\rm Min.}$')
result.rename_column('MAX_M', r'$M_{\rm Max.}$')
result.rename_column('ZMIN', r'$z_{\rm Min.}$')
result.rename_column('ZMAX', r'$z_{\rm Max.}$')
result.rename_column('ZLIMS_NGAL', r'$N_{GAL} / 10^3$')
result.rename_column('VZ', r'$V_{\rm DDP}$ / 10^6$')
result.rename_column('DENS', r'$\rho_{\rm DDP} / 10^{-3}$')

result.pprint()

# https://arxiv.org/pdf/1409.4681.pdf
ascii.write(result, 'tables/Tab2.tex', Writer=ascii.Latex, latexdict=ascii.latex.latexdicts['AA'], overwrite=True)

exit(0)

rows = []

rpath = findfile(ftype='randoms_bd_ddp_n8', survey=survey, version=version, field='G9')
rand  = Table.read(rpath)

print('\n\n')

for idx in np.arange(9):    
    nd8 = 0 
    
    for field in ['G9', 'G12', 'G15']:
        
        fpath = findfile(ftype='ddp_n8_d0', survey=survey, version=version, field=field, utier=idx)
        dat   = Table.read(fpath)
                
        nd8  += len(dat) / 1.e3
        
    rows.append(('d{}'.format(idx), d8_limits[idx][0], d8_limits[idx][1], nd8, rand.meta['DDP1_d{}_VOLFRAC'.format(idx)]))
    
print('\n\n')

# Generate Table 3 of McNaught-Roberts (2014).
result = Table(rows=rows, names=['Label', 'Min_{d8}', 'Max_{d8}', 'N_{d8} [1e3]', 'fd8']) #, 'N_{d8}/N_{max}'])

# TODO: CHECK THE MATHS BY HAND
result['N_{d8} / N_{max}'] = result['N_{d8} [1e3]'] / max(result['N_{d8} [1e3]'])

# BUG: result['TMR N/N_{max}'] = tmr_Nd8 / tmr_Nd8[5]
result['TMR N/N_{max}'] = tmr_Nd8 / tmr_Nd8[4]

for col in result.itercols():
    if col.info.dtype.kind == 'f':        
        np.around(col, decimals=3, out=col)

result.pprint()

opath = 'tables/Tab3.fits'
result.write(opath, format='fits', overwrite=True)

# https://arxiv.org/pdf/1409.4681.pdf
ascii.write(result, 'tables/Tab3.tex', Writer=ascii.Latex, latexdict=ascii.latex.latexdicts['AA'], overwrite=True)

# subprocess.run('cd {}/tables; pdflatex compile_tables.tex'.format(os.environ['CODE_ROOT']), shell=True, check=True)

print('\n\nDone.\n\n')
