import time
import warnings

from   astropy.io.fits.verify import VerifyWarning

# Suppress verify warnings, e.g. HIERARCH card length. 
warnings.simplefilter('ignore', category=VerifyWarning)

def calc_runtime(start, log=None):
    runtime  = time.time() - start
    runtime /= 60.

    if log != None:
        print('{} after {:.4f} mins.'.format(log, runtime))
        
    return  runtime


if __name__ == '__main__':
    start = time.time()

    runtime  = calc_runtime(start, 'Read randoms')
