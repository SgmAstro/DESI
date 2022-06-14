import os
import glob
import astropy.io.fits as     fits

from   pathlib         import Path
from   astropy.io.fits import FITSDiff


def write_filelist(olist, fpath):
  textfile = open(fpath, 'w')

  for ee in olist:
    textfile.write(ee + '\n')
    
  textfile.close()


def run_diff(new_version=4, over=None):
  # https://docs.astropy.org/en/stable/io/fits/api/diff.html 
  if over == None:
      over   = new_version - 1 

  root    = os.environ['TILING_CATDIR']

  new_dir = root + '/v{}/'.format(new_version)
  old_dir = root + '/v{}/'.format(over)

  npaths  = sorted(glob.glob(new_dir + '/*.fits'))
  opaths  = sorted(glob.glob(old_dir + '/*.fits')) 

  print('\n\nDiffing {} files in {} to {} files in {}.\n\n'.format(len(npaths), new_dir, len(opaths), old_dir))

  common  = list(set([os.path.basename(x) for x in npaths]).intersection([os.path.basename(x) for x in opaths]))

  rundiff_dir = '{}/rundiffs/v{}/'.format(os.environ['TILING_CATDIR'], new_version)

  Path(rundiff_dir).mkdir(parents=True, exist_ok=True)

  for cc in common:
    fd     = FITSDiff(new_dir + cc,\
                      old_dir + cc,\
                      ignore_hdus=[],\
                      ignore_keywords=[],\
                      ignore_comments=[],\
                      ignore_fields=[],\
                      numdiffs=10,\
                      rtol=1.e-6,\
                      atol=1.e-6,\
                      ignore_blanks=True,\
                      ignore_blank_cards=True)

    print('{}/{}.txt'.format(rundiff_dir, cc.split('.')[0]))

    fd    = fd.report('{}/{}.txt'.format(rundiff_dir, cc.split('.')[0]), overwrite=True)

  new     = [x for x in npaths if os.path.basename(x) not in common]
  missing = [x for x in opaths if os.path.basename(x) not in common]

  write_filelist(new, '{}/new.txt'.format(rundiff_dir))
  write_filelist(missing, '{}/missing.txt'.format(rundiff_dir))

  print('\n\nDone.\n\n')

if __name__ == '__main__':
    fds = run_diff()
