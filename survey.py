def survey_specifics(survey):
    gama      = {'rlim': 19.8, 'rmax': 12}
    desi      = {'rlim': 19.5, 'rmax': 12}
    
    specifics = {'gama': gama, 'desi': desi}

    assert  survey in specifics.keys()

    return  specifics[survey]
