def survey_specifics(survey):
    gama      = {'rlim': 19.8, 'rmax': 12}

    specifics = {'gama': gama}

    assert  survey in specifics.keys()

    return  specifics[survey]
