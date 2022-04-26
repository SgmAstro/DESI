def survey_specifics(survey):
    
    gama      = {'rlim': 19.8, 'rmax': 12, 'area':         180.}
    desi      = {'rlim': 19.5, 'rmax': 12, 'area': 37.7424 * 6., 'pet_offset': 0.12}
    
    specifics = {'gama': gama, 'desi': desi}

    
    assert  survey in specifics.keys(), f'Requested {survey} is not available.'

    return  specifics[survey]
