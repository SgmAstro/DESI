import numpy as np

# from https://github.com/michaelJwilson/lumfn/blob/master/py/lumfn/GAMA4/in_gama.py
gama_limits = {'G9':  {'ra_min': 129,   'ra_max': 141,   'dec_min': -2, 'dec_max': 3},
               'G12': {'ra_min': 174,   'ra_max': 186,   'dec_min': -3, 'dec_max': 2},
               'G15': {'ra_min': 211.5, 'ra_max': 223.5, 'dec_min': -2, 'dec_max': 3}}

def gama_field(ras, decs):
    result = np.array(['None'] * len(ras))

    for field in gama_limits.keys():
        ra_min = gama_limits[field]['ra_min']
        ra_max = gama_limits[field]['ra_max']

        dec_min = gama_limits[field]['dec_min']
        dec_max = gama_limits[field]['dec_max']

        in_field = (ras >= ra_min) & (ras <= ra_max) & (decs >= dec_min) & (decs <= dec_max)

        result[in_field] = field

    return result
