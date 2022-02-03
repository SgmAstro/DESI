import os
import sys

from   findfile import file_check, fields


DRYRUN  = os.environ['DRYRUN']
GOLDDIR = os.environ['GOLD_DIR']
RANDDIR = os.environ['RANDOMS_DIR']

def tidyup():
    # File check summary. 
    sys.stdout = open(GOLDDIR + 'summary.log', 'w')

    file_check()

    sys.stdout.close()

    # Gather Randoms and write to disk. 
    #
    # fpaths   = findfile('ddp_n8_d0', dryrun=False, prefix='', field=fields, utier=6)
    # all_cats = gather_cat(fpaths)
    # all_cats.pprint()

    
if __name__ == '__main__':
    tidyup()

