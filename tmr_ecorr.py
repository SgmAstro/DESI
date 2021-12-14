import numpy as np

from data.ke_params import *

def tmr_ecorr(redshift, restframe_colour, all=False):
    if all:
        result = -Qall*redshift
    
    else:
        result = np.zeros_like(redshift)
        result[restframe_colour > redblue_split] = -Qred*redshift[restframe_colour > redblue_split]
        result[restframe_colour <= redblue_split] = -Qblue*redshift[restframe_colour <= redblue_split]

    return result

def tmr_q(redshift, restframe_colour, all=False):
    if all:
        result = np.ones_like(redshift) * Qall

    else:
        result = np.zeros_like(redshift)
        result[restframe_colour > redblue_split] = Qred
        result[restframe_colour <= redblue_split] = Qblue

    return result


if __name__ == '__main__':
    import numpy as np

    
    zs = np.arange(0.0, 1.0, 0.01)
    cols = 0.5 * np.ones_like(zs)

    Es = tmr_ecorr(zs, cols, all=False)
    Qs = tmr_q(zs, cols, all=True) 
    
    print(Qs)
