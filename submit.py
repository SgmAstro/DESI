import os
import argparse
import numpy as np


# /cosma/home/durham/dc-wils7/DESI/bin/gold_pipeline
parser = argparse.ArgumentParser(description='Customise pipeline submission scripts.')
parser.add_argument('-s', '--script',  help='Script to customise.',    type=str, default='gold_pipeline')
parser.add_argument('-q', '--queue',   help='Queue for submission.',   type=str, default='cosma')
parser.add_argument('-m', '--memory',  help='Node memory usage [GB].', type=str, default='50')
parser.add_argument('-t', '--time',    help='Job time to request.',    type=str, default='04:00:00')
parser.add_argument('-l', '--log',     help='Job log path.',           type=str, default=None)
parser.add_argument('-a', '--account', help='Account for submission.', type=str, default='desi')
parser.add_argument('-n', '--nodes',   help='Nodes to request.',       type=int, default=2)

args     = parser.parse_args()
script   = args.script
queue    = args.queue
memory   = args.memory
time     = args.time
log      = args.log
account  = args.account
nodes    = args.nodes

if log == None:
    #  findfile
    log  = os.environ['HOME'] + f'/data/GAMA4/logs/{script}.blah'
    
#
ff       = open(f'/cosma/home/durham/dc-wils7/DESI/bin/{script}')
ff       = ff.read()
ff       = ff.split('\n')
ff       = [x.rstrip() for x in ff]

custom   = [x for x in ff if '#SBATCH' in x]
rest     = [x for x in ff if '#SBATCH' not in x] 

#SBATCH -p cordelia
#SBATCH --mem=20G
#SBATCH -t 02:00:00
#SBATCH -o /cosma/home/durham/dc-wils7/data/GAMA4/logs/gold_pipeline.log
#SBATCH -A durham
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --open-mode=append

def swap(xx):
    args = {'queue':   {'sbatch': ['-p', '--partition'], 'arg': queue,        'equals': False},
            'memory':  {'sbatch': ['--mem'],             'arg': memory + 'G', 'equals': True},
            'time':    {'sbatch': ['-t'],                'arg': time,         'equals': False},
            'account': {'sbatch': ['-A'],                'arg': account,      'equals': False},
            'log':     {'sbatch': ['-o'],                'arg': log,          'equals': False},
            'nodes':   {'sbatch': ['--nodes'],           'arg': nodes,        'equals': True}}

    print('\n\n...  Swapping  ...')

    for i, _ in enumerate(xx):
        for arg in sorted(args.keys()):
            var    = args[arg]['arg']
            batch  = np.atleast_1d(args[arg]['sbatch'])
            equals = args[arg]['equals'] 

            if (var != None):
                split = []

                for yy in _.split():
                    split += yy.split('=')

                found = np.any([xx in batch for xx in split])

                # print(var, batch, found, split)

                if found:
                    if equals:
                        xx[i] = '#SBATCH {}={}'.format(batch[0], var)

                    else:
                        xx[i] = '#SBATCH {} {}'.format(batch[0], var)

    return  xx

print('\n\n----  INPUT  ----\n')

for xx in custom:
    print(xx)

##  Swap
custom = swap(custom)

print('\n\n----  OUTPUT  ----\n')

for xx in custom:
    print(xx)

with open(f'{script}.log.tmp', 'w') as f:
    rest.remove('#!/bin/bash')

    to_write = custom + rest

    f.write('#!/bin/bash')
    f.write('\n')

    for line in custom + rest:
        f.write(line)
        f.write('\n')

print('\n\nDone.\n\n')
