import os
import time
import psutil
import warnings

from   astropy.io.fits.verify import VerifyWarning


# Suppress verify warnings, e.g. HIERARCH card length. 
warnings.simplefilter('ignore', category=VerifyWarning)

def calc_runtime(start, log=None, memuse=True):
    runtime  = time.time() - start
    runtime /= 60.

    if memuse:
        process = psutil.Process(os.getpid())
        memuse  = process.memory_info().rss / 1.e6 # [MB]
        memuse  = ' (memuse: {:.2f}MB)'.format(memuse)

    else:
        memuse  = ''

    if log != None:
        print('{} after {:.4f} mins{}.'.format(log, runtime, memuse))
        
    return  runtime


if __name__ == '__main__':
    start = time.time()

    runtime  = calc_runtime(start, 'Read randoms')
