import os


def findfile(ftype, dryrun=False, prefix='', field=None)
    if dryrun:
        dryrun = '_dryrun'
    else:
        dryrun = ''

    gold_dir   = os.environ['GOLD_DIR']
    rand_dir   = os.environ['RANDOMS_DIR']
    
    if field is not None:
        file_types = {'gold':   {'dir': gold_dir, 'id': 'gama_gold', 'ftype': ''},\
                      'kE':     {'dir': gold_dir, 'id': 'gama_gold', 'ftype': 'kE'},\
                      'zmax':   {'dir': gold_dir, 'id': 'gama_gold', 'ftype': 'zmax'},\ 
                      'vmax':   {'dir': gold_dir, 'id': 'gama_gold', 'ftype': 'vmax'},\
                      'lumfn':  {'dir': gold_dir, 'id': 'gama_gold', 'ftype': 'lumfn'},\
                      'ddp':    {'dir': gold_dir, 'id': 'gama_gold', 'ftype': 'ddp'},\
                      'ddp_n8': {'dir': gold_dir, 'id': 'gama_gold', 'ftype': 'ddp_n8'}}

        parts      = file_types[ftype]
        fpath      = parts['dir'] + '{}_{}{}.fits'.format(parts['id'], parts['ftype'], dryrun)

    else:
        file_types = {'ddp_n8_d0': {'dir': gold_dir, 'id': 'gama_gold', 'ftype': ''}

        fpath      = parts['dir'] + '{}_{}_{}{}.fits'.format(parts['id'], field, parts['ftype'], dryrun)
        
    return  fpath
