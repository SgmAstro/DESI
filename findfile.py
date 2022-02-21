import os
import time
import glob
import datetime
import numpy as np

from   astropy.table import Table, vstack
from   delta8_limits import d8_limits
from   gama_limits   import gama_fields
from   desi_fields   import desi_fields

supported = ['gold',\
             'kE',\
             'zmax',\
             'vmax',\
             'lumfn',\
             'ddp',\
             'ddp_n8']

def gather_cat(fpaths):
    if len(fpaths) == 0:
        return  None

    assert  np.all(np.array([os.path.isfile(x) for x in fpaths]))

    tables      = [Table.read(x) for x in fpaths]
    tables      = vstack(tables)

    # TODO:  Headers, e.g. Area, ngal etc.  
    tables.meta = {}

    return  tables 

def fetch_fields(survey):
    if survey == 'gama':
        fields = gama_fields   

    elif survey == 'desi':
        fields = desi_fields

    else:
        raise NotImplementedError

    return fields

def release_dir(user=os.environ['USER'], versioned=True, survey='gama'):
    assert survey == 'gama', 'TODO: Support DESI.'

    # E.g.  /cosma/home/durham/dc-wils7/data/GAMA4/                                                                                                                                                
    if versioned:
        ff = glob.glob('/cosma/home/durham/{}/data/v*'.format(user))
        ff.sort(key=os.path.getmtime)

        return  ff[-1]

    else:
        return '/cosma/home/durham/{}/data/GAMA4/'.format(user)

def overwrite_check(opath):
    if os.path.isfile(opath):
        print('{} found on disk and overwrite forbidden (--nooverwrite).'.format(opath))
        exit(0)
        
def findfile(ftype, dryrun=False, prefix=None, field=None, utier='{utier}', survey='gama', realz=0, debug=True):    
    survey = survey.lower()

    # Special case:                                                                                                                                                                                 
    if (ftype == 'gold') & dryrun & (survey == 'gama'):
        return  os.environ['CODE_ROOT'] + '/data/gama_gold_dryrun.fits'

    if survey == 'gama':
        fields = gama_fields   

    elif survey == 'desi':
        fields = desi_fields

    else:
        raise NotImplementedError
    
    if dryrun:
        dryrun = '_dryrun'

    else:
        dryrun = ''

    realz      = str(realz)

    gold_dir   = os.environ['GOLD_DIR']
    rand_dir   = os.environ['RANDOMS_DIR']

    if isinstance(field, list):
        return  [findfile(ftype, dryrun=dryrun, prefix=prefix, field=ff, utier=utier) for ff in field]
        
    if field == None:
        file_types = {'gold':   {'dir': gold_dir, 'id': f'{survey}',      'ftype': 'gold'},\
                      'kE':     {'dir': gold_dir, 'id': f'{survey}_gold', 'ftype': 'kE'},\
                      'zmax':   {'dir': gold_dir, 'id': f'{survey}_gold', 'ftype': 'zmax'},\
                      'vmax':   {'dir': gold_dir, 'id': f'{survey}_gold', 'ftype': 'vmax'},\
                      'lumfn':  {'dir': gold_dir, 'id': f'{survey}_gold', 'ftype': 'lumfn'},\
                      'ddp':    {'dir': gold_dir, 'id': f'{survey}_gold', 'ftype': 'ddp'},\
                      'ddp_n8': {'dir': gold_dir, 'id': f'{survey}_gold', 'ftype': 'ddp_n8'}}

        parts      = file_types[ftype]
        fpath      = parts['dir'] + '/{}_{}{}.fits'.format(parts['id'], parts['ftype'], dryrun)

    else: 
        file_types = {'ddp_n8_d0':          {'dir': gold_dir, 'id': f'{survey}_gold',         'ftype': 'ddp_n8_d0_{}'.format(utier)},\
                      'ddp_n8_d0_vmax':     {'dir': gold_dir, 'id': f'{survey}_gold',         'ftype': 'ddp_n8_d0_{}_vmax'.format(utier)},\
                      'ddp_n8_d0_lumfn':    {'dir': gold_dir, 'id': f'{survey}_gold',         'ftype': 'ddp_n8_d0_{}_lumfn'.format(utier)},\
                      'randoms':            {'dir': rand_dir, 'id': 'randoms',                'ftype': realz},\
                      'randoms_n8':         {'dir': rand_dir, 'id': 'randoms_N8',             'ftype': realz},\
                      'randoms_bd':         {'dir': rand_dir, 'id': 'randoms_bd',             'ftype': realz},\
                      #'randoms_ddp1':       {'dir': rand_dir, 'id': 'randoms_ddp1',           'ftype': realz},\
                      #'randoms_ddp1_n8':    {'dir': rand_dir, 'id': 'randoms_ddp1_N8',        'ftype': realz},\
                      #'randoms_ddp1_bd':    {'dir': rand_dir, 'id': 'randoms_ddp1_bd',        'ftype': realz},\
                      #'randoms_ddp1_bd_n8': {'dir': rand_dir, 'id': 'randoms_ddp1_bd_ddp_n8', 'ftype': realz},\
                      'randoms_bd_ddp_n8':  {'dir': rand_dir, 'id': 'randoms_bd_ddp_n8',      'ftype': realz}
                     }
        
        parts      = file_types[ftype]
        fpath      = f'' + parts['dir'] + '/{}_{}_{}{}.fits'.format(parts['id'], field, parts['ftype'], dryrun)

        if prefix != None:
            assert 'randoms' in prefix;
            
            dirname = os.path.dirname(fpath)
            fpath = os.path.basename(fpath)
            
            fpath = fpath.replace('randoms', prefix)
            fpath = dirname + '/' + fpath
        
        
    if debug:
        print(f'DEBUG: FINDFILE RETURNS {fpath}')
        
    return  fpath

def file_check(dryrun=None):
    try:
        dryrun = os.environ['DRYRUN']

    except Exception as E:
        print(E)

        dryrun = ''

    fpaths = [findfile(xx, dryrun=False, prefix='', field=None) for xx in supported]

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

    # print('\n\n----  Multiple fields    ----\n')

    return  ~np.all([os.path.isfile(fp) for fp in all_paths])


if __name__ == '__main__':
    failure = file_check()

    print('\n\nSuccess: {}\n\n'.format(~failure))
