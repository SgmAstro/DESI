import sys
import yaml
import inflect

from   collections import OrderedDict


class Configuration:
    name = 'Configuration' # class attributes

    def __init__(self, fpath):
        self.fpath = fpath

        # Load the config data into memory.
        with open(fpath) as f:
            config = yaml.safe_load(f)

        self.config_types = config.keys() 
        self.config_types_attributes = OrderedDict()
        
        for key in self.config_types:
            sub_config = config[key]
            keys       = sub_config.keys()

            self.config_types_attributes[key] = keys

            if key == 'comments':
                self.comments = OrderedDict() # instance attributes

                ckeys  = sub_config.keys() 
                
                for key in ckeys:
                    self.comments[key] = sub_config[key]

            else:
                for key in keys:
                    setattr(self, key, sub_config[key])

    def update_comments(self, comments):
        p  = inflect.engine()
        ps = [p.ordinal(i) for i in range(1, 50, 1)] 
        ps = [str(p) for p in ps if p not in self.comments.keys()]

        for p, comment in zip(ps, comments):        
            self.comments[p] = comment
        
    def get_attributes(self):
        for key in self.config_types:
            if key == 'comments':
                continue

            keys = self.config_types_attributes[key]

            for kk in keys:
                print('{}\t\t{}\t\t{}'.format(key.ljust(20), kk.ljust(20),  getattr(self, kk)))

            

    def write(self, opath):
        keep       = sys.stdout

        sys.stdout = open(opath, 'w')

        keys = config.comments.keys()

        print('Comments:\n')
        
        for kk in keys:
            print(kk, config.comments[kk])

        print('\n\n')

        config.get_attributes()

        sys.stdout.close()
        
        sys.stdout = keep

if __name__ == '__main__':
    config = Configuration('configs/config.yaml')
    config.update_comments(['New comment'])
    
    config.write('/cosma/home/durham/dc-wils7/DESI/configs/config.txt')

    print('\n\nDone.\n\n')
