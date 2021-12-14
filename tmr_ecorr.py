from data.ke_params import *

def tmr_ecorr(redshift, restframe_colour, all=False):
    if all:
        result = -Qall*redshift
    
    else:
        result = np.zeros_like(redshift)
        result[restframe_colour > redblue_split] = -Qred*redshift[restframe_colour > redblue_split]
        result[restframe_colour <= redblue_split] = -Qblue*redshift[restframe_colour <= redblue_split]

    return result

if __name__ == '__main__':
    import numpy as np
    
    zs = np.arange(0.0, 1.0, 0.01)
    cols = 0.5 * np.ones_like(zs)

    Es = tmr_ecorr(zs, cols, all=False)

    print(Es)
