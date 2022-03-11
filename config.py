class config:
    def __init__(self, ff=False, self_count=False):
        self.vmax_fillfactor = ff
        self.self_count = self_count
        
    def vmax(self):
        return self.vmax_fillfactor
        
    def sc(self):
        return self.self_count
    

configuration = config(ff=True, self_count=True)

configuration.sc()