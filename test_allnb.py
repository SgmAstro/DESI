import argparse
import papermill as pm

from   tidyup   import   tidyup
from   findfile import   fetch_fields

parser  = argparse.ArgumentParser(description='Select field.')
parser.add_argument('-s', '--survey', help='Survey, e.g. GAMA, DESI, etc.', type=str, default='gama')

args    = parser.parse_args()
survey  = args.survey.lower()

fields = fetch_fields(survey)

# https://docs.pytest.org/en/6.2.x/
def test_allnbs(survey=survey):
    print('Running all tests.')
    
    survey = survey.lower()
    
    # TODO: Check this
    if survey != 'gama' and survey != 'desi':
        raise NotImplementedError(f'No implementation for survey: {survey}')
    
    tidyup()
    
    run_randomqa(survey)

    run_goldqa(survey)

    run_delta8qa(survey)
    
    print('Done.')
    
def run_randomqa(survey):
    for field in fields:
        print('Running random QA for field {}'.format(field))
        
        pm.execute_notebook('docs/nb/randoms_n8_qa.ipynb',\
                            'test/pm_randoms_n8_{}_qa.ipynb'.format(field),\
                            parameters=dict(field=field,survey=survey),\
                            kernel='lumfn',\
                            )

def run_goldqa(survey):
    print('Running gold QA')

    tests = ['zmax_catQA', 'kE_catQA', 'ddp_QA', 'lumfn', 'delta8_qa']

    for test in tests:
        try:
            pm.execute_notebook('docs/nb/{}.ipynb'.format(test),\
                                'test/pm_{}.ipynb'.format(test),\
                                parameters=dict(survey=survey),\
                                kernel='lumfn',\
            )

        except Error as E:
            print('Failed on {} test with error: '.format(test))
            print(E)

    '''
    pm.execute_notebook('docs/nb/kE_catQA.ipynb',\
                        'test/pm_kE_catQA.ipynb',\
                        kernel='lumfn',\
                        )

    pm.execute_notebook('docs/nb/ddp_QA.ipynb',\
                        'test/pm_ddp_QA.ipynb',\
                        kernel='lumfn',\
                        )

    pm.execute_notebook('docs/nb/lumfn.ipynb',\
                        'test/pm_lumfn.ipynb',\
                        kernel='lumfn',\
                        )

    pm.execute_notebook('docs/nb/delta8_qa.ipynb',\
                        'test/pm_delta8_qa.ipynb',\
                        kernel='lumfn',\
                        )
    '''
    # jack knife qa. 
    # d8 LF. 
    # desi qa. 
    
def run_delta8qa(survey):
    for field in fields:
        print('Running delta8 QA for field {}'.format(field))

        pm.execute_notebook('docs/nb/delta8_qa.ipynb',\
                            'test/pm_delta8_qa_{}.ipynb'.format(field),\
                            parameters=dict(field=field, survey=survey),\
                            kernel='lumfn',\
                            )