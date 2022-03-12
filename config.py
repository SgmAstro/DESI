import os
import sys
import yaml
import inflect
import contextlib

from   pathlib import Path
from   collections import OrderedDict


class CustomDumper(yaml.SafeDumper):
    # inspired by https://stackoverflow.com/a/44284819/3786245
    def write_line_break(self, data=None):
        super().write_line_break(data)

        if len(self.indents) == 1:
            super().write_line_break()

@contextlib.contextmanager
def smart_open(filename=None):
    if filename and filename != '-':
        fh = open(filename, 'w')
    else:
        fh = sys.stdout

    try:
        yield fh
    finally:
        if fh is not sys.stdout:
            fh.close()

class Configuration:
    name = 'Configuration'

    def __init__(self, fpath):
        if fpath == None:
            fpath = os.environ['CODE_ROOT'] + '/configs/config.yaml'

        self.fpath = fpath

        print(f'Reading {fpath}')

        # Load the config data into memory.
        with open(fpath) as f:
            config = yaml.safe_load(f)

        self.config_types = sorted(list(config.keys())) 
        self.config_types.remove('runtime')
        self.config_types.remove('comments')
        
        self.config_types = ['runtime'] + self.config_types + ['comments']
        self.attributes   = OrderedDict()

        # print(self.config_types)
        
        for key in self.config_types:
            self.attributes[key] = config[key]

            sub_config = config[key]
            keys       = sub_config.keys()

            for key in keys:
                setattr(self, key, sub_config[key])

        if self.attributes['runtime']['replay']:
            raise NotImplementedError()

    def update_comments(self, comments):
        p  = inflect.engine()
        ps = [p.ordinal(i) for i in range(1, 50, 1)] 
        ps = [str(p) for p in ps if p not in self.attributes['comments'].keys()]

        for p, comment in zip(ps, comments):        
            self.attributes['comments'][p] = comment

    def print_attributes(self, output=None):
        attr = dict(self.attributes)

        # https://stackoverflow.com/questions/17602878/how-to-handle-both-with-open-and-sys-stdout-nicely
        with smart_open(output) as fh:
            yaml.dump(attr, fh, default_flow_style=False, Dumper=CustomDumper, sort_keys=False)

    def update_attributes(self, new):
        if hasattr(new, '__dict__'):
            new = new.__dict__

        assert isinstance(new, dict)

        for key in sorted(new.keys()):
            setattr(self, key, new[key])
        
        for key in self.config_types:
            keys = self.attributes[key].keys()

            for subkey in keys:
                self.attributes[key][subkey] = getattr(self, subkey)

    def setup_replay(self, new_vars):
        for key in self.config_types:
            new_vars.update(self.attributes[key])

    def write(self, opath):
        self.attributes['runtime']['user'] = os.environ['USER']

        Path(os.path.dirname(opath)).mkdir(parents=True, exist_ok=True)

        print(f'Writing {opath}') 
        
        with open(opath, 'w', encoding = 'utf-8') as ofile:
            yaml.dump(self.attributes, ofile, default_flow_style=False, Dumper=CustomDumper)


if __name__ == '__main__':
    config = Configuration('configs/config.yaml')
    config.update_comments(['New comment'])

    # TODO findfile.
    opath  = os.environ['GOLD_DIR'] + '/configs/config.yaml'
    
    config.write(opath)

    print('\n\nDone.\n\n')
