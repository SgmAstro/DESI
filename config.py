import yaml

class config:
    def __init__(self, fpath):
        self.fpath = fpath

        # Load the config data into memory.
        with open(fpath) as f:
            config = yaml.safe_load(f)

        self.config_types = config.keys() 
        
        for key in self.config_types:
            sub_config = config[key]

            keys = sub_config.keys()

            for key in keys:
                setattr(self, key, sub_config[key])

    def write(self, opath):
        '''
        with open(opath, 'w') as file:
            documents = yaml.dump(, file)
        '''

configuration = config('configs/config.yaml')
