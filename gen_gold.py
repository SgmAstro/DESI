import argparse

from   gama_gold import gama_gold

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Gen kE cat.')
    parser.add_argument('--survey',       help='Survey', default='gama')
    parser.add_argument('--nooverwrite',  help='Do not overwrite outputs if on disk', action='store_true')

    args   = parser.parse_args()

    survey = args.survey

    if survey == 'gama':
        gama_gold(args)

    elif survey == 'desi':
        raise NotImplementedError()
    
    else:
        raise ValueError(f'Survey: {survey} is not supported.')
