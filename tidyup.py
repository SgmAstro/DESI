import os
import sys

from   findfile import file_check, findfile

def tidyup():
    # File check summary. 
    fpath      = findfile('summary_log')

    keep       = sys.stdout

    sys.stdout = open(fpath, 'w')

    file_check()

    sys.stdout.close()
    
    sys.stdout = keep

    
if __name__ == '__main__':
    tidyup()

