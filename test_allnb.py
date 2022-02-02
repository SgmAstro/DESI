import papermill as pm

from   tidyup import tidyup


# TODO: Define elsewhere.
fields = ['G9', 'G12', 'G15']

# https://docs.pytest.org/en/6.2.x/
def test_allnbs():
    print('Running all tests.')

    tidyup()
    
    run_randomqa()

    run_goldqa()

    run_delta8qa()
    
    print('Done.')
    
def run_randomqa():
    for field in fields:
        print('Running random QA for field {}'.format(field))
        
        pm.execute_notebook('docs/nb/randoms_n8_qa.ipynb',\
                            'test/pm_randoms_n8_{}_qa.ipynb'.format(field),\
                            parameters=dict(field=field),\
                            kernel='lumfn',\
                            )

def run_goldqa():
    print('Running gold QA')

    tests = ['zmax_catQA', 'kE_catQA', 'ddp_QA', 'lumfn', 'delta8_qa']

    for test in tests:
        try:
            pm.execute_notebook('docs/nb/{}.ipynb'.format(test),\
                                'test/pm_{}.ipynb'.format(test),\
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
    
def run_delta8qa():
    for field in fields:
        print('Running delta8 QA for field {}'.format(field))

        pm.execute_notebook('docs/nb/delta8_qa.ipynb',\
                            'test/pm_delta8_qa_{}.ipynb'.format(field),\
                            parameters=dict(field=field),\
                            kernel='lumfn',\
                            )
