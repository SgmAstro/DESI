def survey_specifics(survey):
    '''
    Notes:
        Max. separation:  expected maximum separation for boundary dist [Mpc/h].
    '''
    gama      = {'rlim': 19.8, 'rmax': 12, 'area':         180., 'max_sep': 70.}
    desi      = {'rlim': 19.5, 'rmax': 12, 'area': 37.7424 * 6., 'pet_offset': 0.12, 'max_sep': 10.}
    
    specifics = {'gama': gama,\
                 'desi': desi}

    assert  survey in specifics.keys(), f'Requested {survey} is not available.'

    return  specifics[survey]
