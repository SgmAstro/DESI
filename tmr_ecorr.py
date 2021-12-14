def tmr_ecorr(redshift, restframe_colour, all=False):
    if all:
        result = -Q_ALL*redshift
    
    else:
        result = np.zeros_like(redshift)
        result[restframe_colour > 0.62] = -Q_RED*redshift[restframe_colour > 0.62]
        result[restframe_colour <= 0.62] = -Q_BLUE*redshift[restframe_colour <= 0.62]

    return result
