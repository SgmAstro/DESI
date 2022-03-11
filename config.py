import yaml
import inflect

from   collections import OrderedDict


class config:
    def __init__(self, fpath):
        self.fpath = fpath

        # Load the config data into memory.
        with open(fpath) as f:
            config = yaml.safe_load(f)

        self.config_types = config.keys() 
        
        for key in self.config_types:
            sub_config = config[key]
            keys       = sub_config.keys()

            if key == 'comments':
                self.comments = OrderedDict()

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

    def write(self, opath):
        '''
        with open(opath, 'w') as file:
            documents = yaml.dump(, file)
        '''

if __name__ == '__main__':
    configuration = config('configs/config.yaml')
    configuration.update_comments(['New comment'])

    keys = configuration.comments.keys()

    for kk in keys:
        print(kk, configuration.comments[kk])

    print('\n\nDone.\n\n')
