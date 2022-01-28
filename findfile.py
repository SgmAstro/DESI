import os
import time
import glob
import datetime

from   astropy.table import Table, vstack
from   delta8_limits import d8_limits


def findfile(ftype, dryrun=False, prefix='', field=None, utier='{utier}'):
    if dryrun:
        dryrun = '_dryrun'
    else:
        dryrun = ''

    gold_dir   = os.environ['GOLD_DIR']
    rand_dir   = os.environ['RANDOMS_DIR']
    
    if field == None:
        file_types = {'gold':   {'dir': gold_dir, 'id': 'gama',      'ftype': 'gold'},\
                      'kE':     {'dir': gold_dir, 'id': 'gama_gold', 'ftype': 'kE'},\
                      'zmax':   {'dir': gold_dir, 'id': 'gama_gold', 'ftype': 'zmax'},\
                      'vmax':   {'dir': gold_dir, 'id': 'gama_gold', 'ftype': 'vmax'},\
                      'lumfn':  {'dir': gold_dir, 'id': 'gama_gold', 'ftype': 'lumfn'},\
                      'ddp':    {'dir': gold_dir, 'id': 'gama_gold', 'ftype': 'ddp'},\
                      'ddp_n8': {'dir': gold_dir, 'id': 'gama_gold', 'ftype': 'ddp_n8'}}

        parts      = file_types[ftype]
        fpath      = parts['dir'] + '{}_{}{}.fits'.format(parts['id'], parts['ftype'], dryrun)

    else: 
        file_types = {'ddp_n8_d0':          {'dir': gold_dir, 'id': 'gama_gold',              'ftype': 'ddp_n8_d0_{}'.format(utier)},\
                      'ddp_n8_d0_vmax':     {'dir': gold_dir, 'id': 'gama_gold',              'ftype': 'ddp_n8_d0_{}_vmax'.format(utier)},\
                      'ddp_n8_d0_lumfn':    {'dir': gold_dir, 'id': 'gama_gold',              'ftype': 'ddp_n8_d0_{}_lumfn'.format(utier)},\
                      'randoms':            {'dir': rand_dir, 'id': 'randoms',                'ftype':'0'},\
                      'randoms_n8':         {'dir': rand_dir, 'id': 'randoms_N8',             'ftype':'0'},\
                      'randoms_bd':         {'dir': rand_dir, 'id': 'randoms_bd',             'ftype':'0'},\
                      'randoms_ddp1':       {'dir': rand_dir, 'id': 'randoms_ddp1',           'ftype':'0'},\
                      'randoms_ddp1_n8':    {'dir': rand_dir, 'id': 'randoms_ddp1_N8',        'ftype':'0'},\
                      'randoms_ddp1_bd':    {'dir': rand_dir, 'id': 'randoms_ddp1_bd',        'ftype':'0'},\
                      'randoms_ddp1_bd_n8': {'dir': rand_dir, 'id': 'randoms_ddp1_bd_ddp_n8', 'ftype':'0'},\
                      'randoms_bd_ddp_n8':  {'dir': rand_dir, 'id': 'randoms_bd_ddp_n8',      'ftype':'0'}
                     }
        
        parts      = file_types[ftype]

        fpath      = parts['dir'] + '{}_{}_{}{}.fits'.format(parts['id'], field, parts['ftype'], dryrun)
        
    return  fpath


if __name__ == '__main__':
    fields = ['G9', 'G12', 'G15']

    fpaths = []
    fpaths.append(findfile('gold',      dryrun=False, prefix='', field=None))
    fpaths.append(findfile('kE',        dryrun=False, prefix='', field=None))
    fpaths.append(findfile('zmax',      dryrun=False, prefix='', field=None))
    fpaths.append(findfile('vmax',      dryrun=False, prefix='', field=None))
    fpaths.append(findfile('lumfn',     dryrun=False, prefix='', field=None))
    fpaths.append(findfile('ddp',       dryrun=False, prefix='', field=None))
    fpaths.append(findfile('ddp_n8',    dryrun=False, prefix='', field=None))

    for field in fields:
        fpaths.append(findfile('randoms',            dryrun=False, prefix='', field=field))
        fpaths.append(findfile('randoms_n8',         dryrun=False, prefix='', field=field))
        fpaths.append(findfile('randoms_bd',         dryrun=False, prefix='', field=field))
        fpaths.append(findfile('randoms_ddp1',       dryrun=False, prefix='', field=field))
        fpaths.append(findfile('randoms_ddp1_n8',    dryrun=False, prefix='', field=field))
        fpaths.append(findfile('randoms_ddp1_bd',    dryrun=False, prefix='', field=field))
        fpaths.append(findfile('randoms_ddp1_bd_n8', dryrun=False, prefix='', field=field))
        fpaths.append(findfile('randoms_bd_ddp_n8',  dryrun=False, prefix='', field=field))

        for ii, _ in enumerate(d8_limits):
            fpaths.append(findfile('ddp_n8_d0',       dryrun=False, prefix='', field=field, utier=ii))
            fpaths.append(findfile('ddp_n8_d0_vmax',  dryrun=False, prefix='', field=field, utier=ii))
            fpaths.append(findfile('ddp_n8_d0_lumfn', dryrun=False, prefix='', field=field, utier=ii))
    
    
    print('\n\n----  SUPPORTED FPATHS    ----\n')

    for fp in fpaths:
        if os.path.isfile(fp):
            mtime = os.path.getmtime(fp)
            mtime = datetime.datetime.utcfromtimestamp(mtime).strftime('%Y-%m-%d %H:%M:%S')
            
        else:
            mtime = ''

        print('{}\t\t{}\t{}'.format(fp.ljust(100), os.path.isfile(fp), mtime))

    gold_paths  = sorted(glob.glob(os.environ['GOLD_DIR']    + '/*.fits'))
    rand_paths  = sorted(glob.glob(os.environ['RANDOMS_DIR'] + '/*.fits'))
    all_paths   = gold_paths + rand_paths
    
    unsupported = [x for x in all_paths if (x not in fpaths and 'dryrun' not in x)]

    print('\n\n----  UNSUPPORTED FPATHS    ----\n')
    
    for fp in unsupported:
        if os.path.isfile(fp):
            mtime = os.path.getmtime(fp)
            mtime = datetime.datetime.utcfromtimestamp(mtime).strftime('%Y-%m-%d %H:%M:%S')
        else:
            mtime = ''

        print('{}\t\t{}\t{}'.format(fp.ljust(100), os.path.isfile(fp), mtime))

    if len(fp) == 0:
        print('WARNING: NO FILES FOUND!')        
    
    print('\n\nDone.\n\n')
