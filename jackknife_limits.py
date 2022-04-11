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



def jk_zones(dat, rand, jk_ra, jk_dec, plot=False):
    '''
    Splits up single GAMA field into jackknife areas based on randoms.
    '''
    
    assert jk_ra > 1, 'jk_ra must be greater than 1'
    assert jk_dec > 1, 'jk_dec must be greater than 1'
    
    rand.sort(['RANDOM_RA', 'RANDOM_DEC'])
    x = rand['RANDOM_RA']

    rand.sort(['RANDOM_DEC', 'RANDOM_RA'])
    y = rand['RANDOM_DEC']

    bin_x = []
    for idx in range(1, jk_ra):
        bin_x.append(np.percentile(rand['RANDOM_RA'], idx/jk_ra*100))

    bin_y = []
    for idx in range(1, jk_dec):
        bin_y.append(np.percentile(rand['RANDOM_DEC'], idx/jk_dec*100))

    bins = [bin_x, bin_y]
    
    dat['JK_RA'] = np.digitize(dat['RA'],bin_x,right=True)
    dat['JK_DEC'] = np.digitize(dat['DEC'],bin_y,right=True)
    dat['JK'] = dat['JK_RA'] * jk_ra + dat['JK_DEC']
    
    del dat['JK_RA']
    del dat['JK_DEC']
    
    if plot:
        for idx in np.unique(dat['JK']):
            plt.scatter(dat[dat['JK'] == idx]['RA'], dat[dat['JK'] == idx]['DEC'], s=0.25)

        plt.xlabel('RA')
        plt.ylabel('DEC')
        plt.gca().set_aspect("equal")
        plt.show()

    return dat