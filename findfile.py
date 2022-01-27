import os
import glob

from   astropy.table import Table, vstack


def findfile(ftype, dryrun=False, prefix='', field=None):
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
        file_types = {'ddp_n8_d0': {'dir': gold_dir, 'id': 'gama_gold', 'ftype': 'ddp_n8_d0'}}
        parts      = file_types[ftype]

        fpath      = parts['dir'] + '{}_{}_{}{}.fits'.format(parts['id'], field, parts['ftype'], dryrun)
        
    return  fpath


if __name__ == '__main__':

    fpaths = []
    fpaths.append(findfile('gold',   dryrun=False, prefix='', field=None))
    fpaths.append(findfile('kE',     dryrun=False, prefix='', field=None))
    fpaths.append(findfile('zmax',   dryrun=False, prefix='', field=None))
    fpaths.append(findfile('vmax',   dryrun=False, prefix='', field=None))
    fpaths.append(findfile('lumfn',  dryrun=False, prefix='', field=None))
    fpaths.append(findfile('ddp',    dryrun=False, prefix='', field=None))
    fpaths.append(findfile('ddp_n8', dryrun=False, prefix='', field=None))
    
    print('\n\n----  SUPPORTED FPATHS    ----\n')

    for fp in fpaths:
        print('{}\t\t{}'.format(fp.ljust(100), os.path.isfile(fp)))

    gold_paths  = sorted(glob.glob(os.environ['GOLD_DIR']    + '/*.fits'))
    rand_paths  = sorted(glob.glob(os.environ['RANDOMS_DIR'] + '/*.fits'))
    all_paths   = gold_paths + rand_paths
    
    unsupported = [x for x in all_paths if (x not in fpaths and 'dryrun' not in x)]

    print('\n\n----  UNSUPPORTED FPATHS    ----\n')

    for fp in unsupported:
        print('{}\t\t{}'.format(fp.ljust(100), os.path.isfile(fp)))

    print('\n\nDone.\n\n')
