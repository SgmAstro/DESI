import os
import sys
import argparse
import subprocess
import numpy as np

from   astropy.io    import ascii
from   astropy.table import Table
from   ddp           import tmr_DDP1, tmr_DDP2, tmr_DDP3
from   delta8_limits import delta8_tier, d8_limits
from   findfile      import findfile, fetch_header


parser = argparse.ArgumentParser(description='Generate Summary Stats')
parser.add_argument('-s', '--survey', help='Select survey', default='gama')

args             = parser.parse_args()
survey           = args.survey

fpath            = findfile(ftype='ddp', survey=survey)
dat              = Table.read(fpath)
names            = ['ZMIN', 'ZMAX', 'VZ', 'DENS']

tmr_DDPs         = [tmr_DDP1, tmr_DDP2, tmr_DDP3]

result = Table()
rows   = []

print('\n\n')

rpath  = findfile(ftype='randoms_bd_ddp_n8', dryrun=False, field='GALL', survey=survey, prefix='randoms_ddp1')

print(f'Fetching {rpath}')

for ddp, tmr_DDP in zip(np.arange(1, 4, 1), tmr_DDPs):
    name  = 'DDP{}_FULL8FRAC'.format(ddp)
    value = fetch_header(fpath=rpath, name=name)

    row   = [ddp, tmr_DDP[0], tmr_DDP[1]]
    
    for col in names:
        row += [dat.meta['DDP{}_{}'.format(ddp, col)]]
        
    row += [dat.meta['DDP{}{}'.format(ddp, 'ZLIMS_NGAL')]]        
    row += [value]

    row  = tuple(row)         

    rows.append(row)

names  = ['DDP', 'MIN_M', 'MAX_M'] + names + ['ZLIMS_NGAL', 'DDP_FULL8FRAC']
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
result.rename_column('DDP_FULL8FRAC', r'Complete frac.')

result.pprint()

# https://arxiv.org/pdf/1409.4681.pdf
ascii.write(result, 'tables/Tab2.tex', Writer=ascii.Latex, latexdict=ascii.latex.latexdicts['AA'], overwrite=True)

rows  = []

rpath = findfile(ftype='randoms_bd_ddp_n8', dryrun=False, field='GALL', survey=survey, prefix='randoms_ddp1')
rand  = Table.read(rpath)

print('\n\n')

for idx in np.arange(9):    
    nd8 = 0 
    
    for field in ['G9', 'G12', 'G15']:        
        # Number of galaxies in each density tier (across fields). 
        fpath = findfile(ftype='ddp_n8_d0', survey=survey, field=field, utier=idx)
        dat   = Table.read(fpath)
                
        nd8  += len(dat) / 1.e3

    # 3-field. 
    volfrac    = rand.meta['DDP1_d{}_VOLFRAC'.format(idx)] 
    volfrac_zp = rand.meta['DDP1_d{}_ZEROPOINT_VOLFRAC'.format(idx)] 

    d8         = rand.meta['DDP1_d{}_TIERMEDd8'.format(idx)]
        
    rows.append(('d{}'.format(idx), d8_limits[idx][0], d8_limits[idx][1], nd8, volfrac, volfrac_zp, d8))
    
print('\n\n')

# Generate Table 3 of McNaught-Roberts (2014).
names  = np.array(['Label', r'Min. $\delta_8$', r'Max. $\delta_8$', '\# galaxies [1e3]', 'Vol. frac.', 'Zeropoint vol. frac.', r'$\langle \delta_8 \rangle$'])
result = Table(rows=rows, names=names)

for name in names[1:]:
    result[name] = result[name].data.astype(float)
    result[name] = np.round(result[name], 4)

result.pprint()

opath = 'tables/Tab3.fits'
result.write(opath, format='fits', overwrite=True)

# https://arxiv.org/pdf/1409.4681.pdf
ascii.write(result, 'tables/Tab3.tex', Writer=ascii.Latex, latexdict=ascii.latex.latexdicts['AA'], overwrite=True)

print('\n\nDone.\n\n')
