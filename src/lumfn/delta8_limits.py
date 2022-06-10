import numpy as np

# Table 3 of https://arxiv.org/pdf/1409.4681.pdf
d8_limits = [[-1.00, -0.75],\
             [-0.75, -0.55],\
             [-0.55, -0.40],\
             [-0.40,  0.00],\
             [ 0.00,  0.70],\
             [ 0.70,  1.60],\
             [ 1.60,  2.90],\
             [ 2.90,  4.00],\
             [ 4.00,  1.e12]]

d8_limits   = np.array(d8_limits)
d8_plot_idx = [0, 3, 5, 8]

# Derived from our overleaf, approx.
d8_tmr      = [-0.895, -0.657, 0.4986, -0.2156, 0.3287, 1.1053, 2.1505, 3.3827, 5.0958]

def delta8_tier(delta8):
    result = -99 * np.ones(len(delta8), dtype=np.int)

    for i, lims in enumerate(d8_limits):
        result[(delta8 >= lims[0]) & (delta8 < lims[1])] = i

    return  result
