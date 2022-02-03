import os
import time
import glob
import datetime

from   astropy.table import Table, vstack
from   delta8_limits import d8_limits


fields    = ['G9', 'G12', 'G15']
supported = ['gold',\
             'kE',\
             'zmax',\
             'vmax',\
             'lumfn',\
             'ddp',\
             'ddp_n8']

def gather_cat(fpaths):
    tables      = [Table.read(x) for x in fpaths]
    tables      = vstack(tables)

    # TODO:  Headers, e.g. Area, ngal etc.  
    tables.meta = {}

    return  tables 


def findfile(ftype, dryrun=False, prefix='', field=None, utier='{utier}', survey='gama'):
    if dryrun:
        dryrun = '_dryrun'
    else:
        dryrun = ''

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
        fpath      = parts['dir'] + '{}_{}{}.fits'.format(parts['id'], parts['ftype'], dryrun)

    else:
        file_types = {'ddp_n8_d0': {'dir': gold_dir, 'id': f'{survey}_gold', 'ftype': 'ddp_n8_d0_{}'.format(utier)}}
        parts      = file_types[ftype]

        fpath      = f'' + parts['dir'] + '{}_{}_{}{}.fits'.format(parts['id'], field, parts['ftype'], dryrun)

    return  fpath

def file_check(dryrun=None):
    try:
        dryrun = os.environ['DRYRUN']

    except Exception as E:
        print(E)

        dryrun = ''

    fpaths = [findfile(xx, dryrun=False, prefix='', field=None) for xx in supported]

    for field in fields:
        for ii, _ in enumerate(d8_limits):
            fpaths.append(findfile('ddp_n8_d0', dryrun=False, prefix='', field=field, utier=ii))
    
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


    print('\n\n----  Multiple fields    ----\n')

    # return  np.all([os.path.isfile(fp) for fp in fpaths])


if __name__ == '__main__':
    file_check()

    print('\n\nDone.\n\n')
