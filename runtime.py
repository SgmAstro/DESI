import time


def calc_runtime(start, log=None):
    runtime  = time.time() - start
    runtime /= 60.

    if log != None:
        print('{} after {:.4f} mins.'.format(log, runtime))
        
    return  runtime


if __name__ == '__main__':
    start = time.time()

    runtime  = calc_runtime(start, 'Read randoms')
