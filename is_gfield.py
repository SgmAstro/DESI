from gama_limits import gama_limits

def is_gfield(ra, dec, gfield):
    
    if (ra > gama_limits[gfield]['ra_min']) & (ra < gama_limits[gfield]['ra_max']) & (dec > gama_limits[gfield]['dec_min']) & (dec < gama_limits[gfield]['dec_max']):
        return 1
    else:
        return 0