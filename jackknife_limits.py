import numpy as np

jk_limits = {'JK0': {'ra_min': 129., 'ra_max': 133.,   'dec_min': -2., 'dec_max': 3.},
            'JK1': {'ra_min': 133., 'ra_max': 137., 'dec_min': -2., 'dec_max': 3.},
            'JK2': {'ra_min': 137., 'ra_max': 141., 'dec_min': -2., 'dec_max': 3.},
            'JK3': {'ra_min': 174., 'ra_max': 178., 'dec_min': -3., 'dec_max': 2.},
            'JK4': {'ra_min': 178., 'ra_max': 182., 'dec_min': -3., 'dec_max': 2.},
            'JK5': {'ra_min': 182., 'ra_max': 186., 'dec_min': -3., 'dec_max': 2.},
            'JK6': {'ra_min': 211.5, 'ra_max': 215.5, 'dec_min': -2., 'dec_max': 3.},
            'JK7': {'ra_min': 215.5, 'ra_max': 219.5, 'dec_min': -2., 'dec_max': 3.},
            'JK8': {'ra_min': 219.5, 'ra_max': 223.5, 'dec_min': -2., 'dec_max': 3.}
            }


def jk_field(ras, decs):
    result = np.array(['None'] * len(ras), dtype=np.str)

    for strip in jk_limits.keys():
        ra_min = jk_limits[strip]['ra_min']
        ra_max = jk_limits[strip]['ra_max']

        dec_min = jk_limits[strip]['dec_min']
        dec_max = jk_limits[strip]['dec_max']

        in_strip = (ras >= ra_min) & (ras <= ra_max) & (decs >= dec_min) & (decs <= dec_max)

        result[in_strip] = strip

    return result
