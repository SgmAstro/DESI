def vol_normalise(cat, scale):
    cat['PHI_N'] = cat['PHI_N']*scale
    cat['PHI_N_ERROR'] = cat['PHI_N_ERROR']*scale
    cat['PHI_IVMAX'] = cat['PHI_IVMAX']*scale
    cat['PHI_IVMAX_ERROR'] = cat['PHI_IVMAX_ERROR']*scale
    
    return cat
    