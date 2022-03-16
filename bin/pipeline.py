
import  os
import  sys
import  glob
import  argparse
import  subprocess
import  numpy      as     np

from    pathlib    import Path
from    subprocess import check_output
from    findfile   import fetch_fields
from    submit     import customise_script
from    config     import Configuration

def run_command(cmd, noid=False):
    print('Command: {}'.format(cmd))

    cmd = cmd.split()

    env = {}
    # env.update(os.environ)

    # print('Calling ...')

    out = check_output(cmd)
    out = out.decode('utf-8') 
    out = out.replace('\n', '')
    
    # print(out)

    if noid:
        out=0
    
    out = int(out)
    
    return out

# Sbatch: python3 pipeline.py --survey desi --use_sbatch --log
# Head:   python3 pipeline.py --survey desi 
parser  = argparse.ArgumentParser(description='Run Lumfn pipeline')
parser.add_argument('--use_sbatch',   help='Submit via Sbatch', action='store_true')
parser.add_argument('--reset',        help='Reset', action='store_true')
parser.add_argument('--nooverwrite',  help='Do not overwrite outputs if on disk', action='store_true')
parser.add_argument('--dryrun',       help='Dryrun', action='store_true')
parser.add_argument('--survey',       help='Survey', default='gama')
parser.add_argument('--freshclone',   help='Fresh clone', action='store_true')
parser.add_argument('--log',          help='Log stdout.', action='store_true')
parser.add_argument('--custom',       help='Customised submission scripts.', default=True)
parser.add_argument('--config',       help='Path to configuration file', type=str, default=None)
parser.add_argument('--comments',     help='Add comments to README.')

# Customise submission scripts.
parser.add_argument('-s', '--script',  help='Script to customise.',    type=str, default=None)
parser.add_argument('--script_log',    help='Job log path.',           type=str, default=None)
parser.add_argument('-q', '--queue',   help='Queue for submission.',   type=str, default=None)
parser.add_argument('-m', '--memory',  help='Node memory usage [GB].', type=str, default=None)
parser.add_argument('-t', '--time',    help='Job time to request.',    type=str, default=None)
parser.add_argument('-a', '--account', help='Account for submission.', type=str, default=None)
parser.add_argument('-n', '--nodes',   help='Nodes to request.',       type=int, default=None)

args        = parser.parse_args()
use_sbatch  = int(args.use_sbatch)
reset       = args.reset
nooverwrite = args.nooverwrite
dryrun      = args.dryrun
survey      = args.survey
freshclone  = args.freshclone
custom      = args.custom
comments    = args.comments
config      = args.config

if config != None:
    config   = Configuration(config)

if comments != None:
    comments = comments.split(';')
    
    print('Appending comments: {}'.format(comments))

    config.update_comments(comments)

if custom:
    customise_script(args)

    custom = '/custom/'

else:
    custom = ''

if reset & (survey == 'desi'):
    raise  NotImplementedError('Reset not currently supported for DESI.')

if reset & nooverwrite:
    raise  ValueError('No overwrite and reset are incompatible.')

if dryrun:
    dryrun  = '--dryrun' 
else:
    dryrun  = ''
    
if nooverwrite:
    nooverwrite = '--nooverwrite'
else:
    nooverwrite = ''

if args.log:
    sys.stdout  = open('pipeline.log', 'w')

if reset:
    print('\n\n>>>>>  TRASHING GOLD_DIR AND RANDOMS  <<<<<\n\n')

    cmds = []

    for root in [os.environ['GOLD_DIR'], os.environ['RANDOMS_DIR']]:
        cmds.append('rm -f {}/logs/*.log'.format(root))
        # cmds.append('rm -f {}/*_dryrun.fits'.format(root))

    for cmd in cmds:
        print(cmd)

        os.system(cmd)

if reset:
    os.environ['RESET']   = str(1)
else:
    os.environ['RESET']   = str(0)

os.environ['DRYRUN']	  = dryrun
os.environ['SURVEY']      = survey
os.environ['NOOVERWRITE'] = nooverwrite

print('\n\nAssuming $USESBATCH={}'.format(use_sbatch))
print('Assuming $RESET={}'.format(reset))
print('Assuming $SURVEY={}'.format(survey))
print('Assuming $DRYRUN={}'.format(dryrun))
print('Assuming $NOOVERWRITE={}\n\n'.format(nooverwrite))

#  ----  Total of eight jobs, with correct dependency logic  ----                                                                                                                                    
os.environ['RESET'] = '0'

#  ---------------------------------------------
home = os.environ['HOME']

os.chdir(f'{home}')

if freshclone:
   cmds = []

   cmds.append('rm -rf {}/tmp'.format(os.environ['HOME']))
   cmds.append('mkdir -p {}/tmp'.format(os.environ['HOME']))

# HACK!!!!!
   cmds.append('git clone --branch oversampling_testrun https://github.com/SgmAstro/DESI.git {}/tmp/DESI'.format(os.environ['HOME']))
   
   for cmd in cmds:    
       out = run_command(cmd, noid=True)

   code_root = os.environ['CODE_ROOT'] = '{}/tmp/DESI/'.format(os.environ['HOME'])

   os.environ['PATH'] = f':{code_root}/bin/{custom}:' + os.environ['PATH']
   os.environ['PYTHONPATH'] = f'{code_root}/:' + os.environ['PYTHONPATH']

else:
    code_root = '{}/DESI/'.format(os.environ['HOME'])

    os.environ['PATH'] = code_root + f'/bin/{custom}:' + os.environ['PATH']
    os.environ['PYTHONPATH'] = code_root + '/:' + os.environ['PYTHONPATH']

Path(os.environ['GOLD_DIR'] + '/logs/').mkdir(parents=True, exist_ok=True)
Path(os.environ['RANDOMS_DIR'] + '/logs/').mkdir(parents=True, exist_ok=True)

#  ---------------------------------------------
# Generate all steps up to reference LF. 
cmd = 'serialorparallel -p {:d} -e DRYRUN={},RESET={:d},NOOVERWRITE={},SURVEY={} -s gold_pipeline -c {}'.format(int(use_sbatch), dryrun, int(reset), nooverwrite, survey, code_root)

gold_jobid = run_command(cmd)

print('\n>>>>> GOLD JOB ID <<<<<')
print(gold_jobid)
print('\n\n')

#
# https://slurm.schedmd.com/sbatch.html
#

fields = fetch_fields(survey=survey)

rand_jobids     = {}
rand_ddp_jobids = {}

for field in fields:
    # No dependency.  Generate all steps up to random fill factor and bound_dist.      
    # RAND_G9_JOBID=$(serialorparallel  -p $USESBATCH -e FIELD=G9,DRYRUN=$DRYRUN,RESET=$RESET,NOOVERWRITE=$NOOVERWRITE  -s rand_pipeline -c $CODE_ROOT)                                             
    # RAND_G12_JOBID=$(serialorparallel -p $USESBATCH -e FIELD=G12,DRYRUN=$DRYRUN,RESET=$RESET,NOOVERWRITE=$NOOVERWRITE -s rand_pipeline -c $CODE_ROOT)                                             
    # RAND_G15_JOBID=$(serialorparallel -p $USESBATCH -e FIELD=G15,DRYRUN=$DRYRUN,RESET=$RESET,NOOVERWRITE=$NOOVERWRITE -s rand_pipeline -c $CODE_ROOT)     
    
    cmd = 'serialorparallel -p {:d} -e FIELD={},DRYRUN={},RESET={:d},NOOVERWRITE={},SURVEY={} -s rand_pipeline -c {}'
    cmd = cmd.format(int(use_sbatch), field, dryrun, int(reset), nooverwrite, survey, code_root)

    rand_jobids[field] = run_command(cmd)
    
for field in fields:
    # Dependency on $GOLD_JOBID (gold ddp cat generated by gold_pipeline). 
    # Generate ddp1 randoms limited to ddp1 z limits - with corresponding fillfactors, bound_dist etc. 
    # RAND_DDP_G9_JOBID=$(serialorparallel   -p $USESBATCH -e FIELD=G9,DRYRUN=$DRYRUN,RESET=$RESET,NOOVERWRITE=$NOOVERWRITE  -d $GOLD_JOBID -s rand_ddp1_pipeline -c $CODE_ROOT)                    
    # RAND_DDP_G12_JOBID=$(serialorparallel  -p $USESBATCH -e FIELD=G12,DRYRUN=$DRYRUN,RESET=$RESET,NOOVERWRITE=$NOOVERWRITE -d $GOLD_JOBID -s rand_ddp1_pipeline -c $CODE_ROOT)                    
    # RAND_DDP_G15_JOBID=$(serialorparallel  -p $USESBATCH -e FIELD=G15,DRYRUN=$DRYRUN,RESET=$RESET,NOOVERWRITE=$NOOVERWRITE -d $GOLD_JOBID -s rand_ddp1_pipeline -c $CODE_ROOT) 

    cmd = 'serialorparallel -p {:d} -e FIELD={},DRYRUN={},RESET={:d},NOOVERWRITE={},SURVEY={} -d {} -s rand_ddp1_pipeline -c {}'
    cmd = cmd.format(int(use_sbatch), field, dryrun, int(reset), nooverwrite, survey, gold_jobid, code_root)
    rand_ddp_jobids[field] = run_command(cmd)

print('\n\n>>>>> RANDOM JOB IDS <<<<<')
print(rand_jobids)
print(rand_ddp_jobids)
print('\n\n')

rand_d8_jobids     = {}
rand_ddp_d8_jobids = {}

for field in fields:
    # Requires ddp cat, & randoms; no reset required. 
    # RAND_D8_G9_JOBID=$(serialorparallel   -p $USESBATCH -e FIELD=G9,DRYRUN=$DRYRUN,NOOVERWRITE=$NOOVERWRITE   -d $GOLD_JOBID,$RAND_G9_JOBID  -s rand_d8_pipeline -c $CODE_ROOT) 
    # RAND_D8_G12_JOBID=$(serialorparallel  -p $USESBATCH -e FIELD=G12,DRYRUN=$DRYRUN,NOOVERWRITE=$NOOVERWRITE  -d $GOLD_JOBID,$RAND_G12_JOBID -s rand_d8_pipeline -c $CODE_ROOT)
    # RAND_D8_G15_JOBID=$(serialorparallel  -p $USESBATCH -e FIELD=G15,DRYRUN=$DRYRUN,NOOVERWRITE=$NOOVERWRITE  -d $GOLD_JOBID,$RAND_G15_JOBID -s rand_d8_pipeline -c $CODE_ROOT) 
    rand_jobid = rand_jobids[field]

    cmd = 'serialorparallel -p {:d} -e FIELD={},DRYRUN={},RESET={:d},NOOVERWRITE={},SURVEY={} -d {},{} -s rand_d8_pipeline -c {}'
    
    cmd = cmd.format(int(use_sbatch), field, dryrun, int(reset), nooverwrite, survey, gold_jobid, rand_jobid, code_root)

    rand_d8_jobids[field] = run_command(cmd)
    


    # RAND_DDP_D8_G9_JOBID=$(serialorparallel   -p $USESBATCH -e FIELD=G9,DRYRUN=$DRYRUN,NOOVERWRITE=$NOOVERWRITE   -d $GOLD_JOBID,$RAND_DDP_G9_JOBID  -s rand_ddp1_d8_pipeline -c $CODE_ROOT)
    # RAND_DDP_D8_G12_JOBID=$(serialorparallel  -p $USESBATCH -e FIELD=G12,DRYRUN=$DRYRUN,NOOVERWRITE=$NOOVERWRITE  -d $GOLD_JOBID,$RAND_DDP_G12_JOBID -s rand_ddp1_d8_pipeline -c $CODE_ROOT)
    # RAND_DDP_D8_G15_JOBID=$(serialorparallel  -p $USESBATCH -e FIELD=G15,DRYRUN=$DRYRUN,NOOVERWRITE=$NOOVERWRITE  -d $GOLD_JOBID,$RAND_DDP_G15_JOBID -s rand_ddp1_d8_pipeline -c $CODE_ROOT) 

    rand_ddp_jobid = rand_ddp_jobids[field]
    #cmd = f'serialorparallel -p $USESBATCH -e FIELD={field},DRYRUN=$DRYRUN,NOOVERWRITE=$NOOVERWRITE,SURVEY=$SURVEY -d {gold_jobid},{rand_ddp_jobid} -s rand_ddp1_pipeline -c $CODE_ROOT'
    
    cmd = 'serialorparallel -p {:d} -e FIELD={},DRYRUN={},NOOVERWRITE={},SURVEY={} -d {},{} -s rand_ddp1_d8_pipeline -c {}'
    cmd = cmd.format(int(use_sbatch), field, dryrun, nooverwrite, survey, gold_jobid, rand_ddp_jobid, code_root)
    rand_ddp_d8_jobids[field] = run_command(cmd)

print('\n\n>>>>> RANDOM D8 JOB IDS <<<<<')
print(rand_d8_jobids)
print(rand_ddp_d8_jobids)
print('\n\n')

# Requires ddp cat. & random fill factor.                                                                                                                                                            
# Note: runs all fields simultaneously.  
dependencies = ','.join(str(rand_ddp_d8_jobids[field]) for field in fields)

# possibly missing RESET=$RESET
# cmd = 'serialorparallel -p $USESBATCH -e DRYRUN=$DRYRUN,NOOVERWRITE=$NOOVERWRITE,SURVEY=$SURVEY -d {dependencies} -s gold_d8_pipeline -c $CODE_ROOT'
cmd = 'serialorparallel -p {:d} -e DRYRUN={},NOOVERWRITE={},SURVEY={} -d {} -s gold_d8_pipeline -c {}'
cmd = cmd.format(int(use_sbatch), dryrun, nooverwrite, survey, dependencies, code_root)
gold_d8_jobid = run_command(cmd)

print('\n\n>>>>>  GOLD D8 JOB IDS  <<<<<')
print(gold_d8_jobid)
print('\n\n>>>>>  DONE.  <<<<<\n\n')

if args.log:
    sys.stdout.close()
