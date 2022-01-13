import papermill as pm

# https://docs.pytest.org/en/6.2.x/
def test_allnbs():
    print('Running all tests.')
    
    run_randomqa()

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
