import papermill as pm

# https://docs.pytest.org/en/6.2.x/
def test_allnbs():
    print('Running all tests.')
    
    run_randomqa()

    # run_delta8qa()
    
    print('Done.')
    
# TODO:  GitHub Action pytest calls fails due to randoms path not
# existing. 
def run_randomqa():
    fields = ['G9', 'G12', 'G15']
    
    for field in fields:
        print('Running random QA for field {}'.format(field))
        
        pm.execute_notebook('docs/nb/randoms_qa.ipynb',\
                            'test/pm_randoms_qa_{}.ipynb'.format(field),\
                            parameters=dict(field=field)
                            )

def run_delta8qa():
    fields = ['G9', 'G12', 'G15']

    for field in fields:
        print('Running delta8 QA for field {}'.format(field))

        pm.execute_notebook('docs/nb/delta8_qa.ipynb',\
                            'test/pm_delta8_qa_{}.ipynb'.format(field),\
                            parameters=dict(field=field)
                            )
