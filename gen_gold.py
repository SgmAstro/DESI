import os
import argparse
import numpy as np

from   gama_gold     import gama_gold
from   desi_gold     import desi_gold
from   findfile      import findfile
from   astropy.table import Table


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Gen kE cat.')
    parser.add_argument('--log',          help='Create a log file of stdout.', action='store_true')
    parser.add_argument('--survey',       help='Survey', default='gama')
    parser.add_argument('--dryrun',       help='Dryrun', action='store_true')
    parser.add_argument('--nooverwrite',  help='Do not overwrite outputs if on disk', action='store_true')
    parser.add_argument('--in_bgsbright', help='Add flag for IN_BGSBRIGHT', action='store_true')

    args        = parser.parse_args()
    log         = args.log
    survey      = args.survey
    dryrun      = args.dryrun
    nooverwrite = args.nooverwrite
    
    if survey == 'gama':
        gama_gold(args)

    elif survey == 'desi':
        if 'NERSC_HOST' in os.environ.keys():
            # Support to run on nersc only.                                                                                                                                                         
            desi_gold()

        else:
            opath = findfile(ftype='gold', dryrun=False, survey='desi')

            print(f'As you are not running at nersc, the output of this script is assumed to be present at {opath}.')

            # DEBUG/PATCH
            print('WARNING: patching desi gold @ {}'.format(opath))

            dat = Table.read(opath)
            dat['ZSURV'] = dat['ZDESI']
            dat['IN_D8LUMFN'] = np.zeros_like(dat['FIELD'], dtype=int)

            ##  BUG/PATCH/  High completeness part of the rosette                                                                                                                                         
            ##  TODO        Deal with at level of IN_D8LUMFN.                                                                                                                                            
            dat               = dat[(dat['ROS_DIST'].data > 0.5) & (dat['ROS_DIST'].data < 1.5)]
            dat.write(opath, format='fits', overwrite=True)

            # if not os.path.exists(findfile(ftype='gold', dryrun=True, survey='desi')):
            print('PATCH/WARNING: as desi dryrun file is not present, copying to {}.'.format(findfile(ftype='gold', dryrun=True, survey='desi')))

            idx = np.arange(len(dat))
            idx = np.random.choice(idx, size=len(idx), replace=False)

            dat = dat[idx]
            dat = dat[:5000]
                
            dat.write(findfile(ftype='gold', dryrun=True, survey='desi'), format='fits', overwrite=True)
    
    else:
        raise  ValueError(f'Survey: {survey} is not supported.')
