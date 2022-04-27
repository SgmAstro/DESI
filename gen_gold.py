import os
import argparse

from   gama_gold import gama_gold
from   desi_gold import desi_gold
from   findfile  import findfile


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

            if not os.path.exists(findfile(ftype='gold', dryrun=True, survey='desi')):
                print('WARNING: as desi dryrun file is not present, copying to {}.'.format(findfile(ftype='gold', dryrun=True, survey='desi')))

                dat = Table.read(opath)
                dat = dat[:5000]
                
                dat.write(findfile(ftype='gold', dryrun=True, survey='desi'), format='fits', overwrite=True)
    
    else:
        raise  ValueError(f'Survey: {survey} is not supported.')
