import os
import sys

def error_stash(exctype, value, traceback, stash_dir=None):
    if stash_dir == None:
        stash_dir  = os.environ['GOLD_DIR'] + '/errors/'

    os.remove(stash_dir)

    if exctype != None:
        traceback.print_tb(f'{stash_dir}_errors.txt')

    else:
        sys.__excepthook__(exctype, value, traceback)

sys.excepthook = error_stash


if __name__ == '__main__':
    raise RuntimeError('Errors test.')

    print('Done.')
