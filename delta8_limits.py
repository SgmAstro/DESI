import numpy as np

dd8_limits = [[-1.0, -0.75], [-0.4, 0.0], [0.7, 1.6], [4.0, 1.e6]]

def delta8_tier(delta8):
    result = -99 * np.ones(len(delta8), dtype=np.int)

    # Gaps in defined d8 coverage??  See TMR.   
    for i, lims in enumerate(dd8_limits):
        result[(delta8 > lims[0]) & (delta8 <= lims[1])] = i

    return  result
