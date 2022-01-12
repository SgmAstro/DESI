import papermill as pm

# https://docs.pytest.org/en/6.2.x/
def test_notebooks():
    fields = ['G9', 'G12', 'G15']
    
    for field in fields:    
        pm.execute_notebook('docs/nb/randoms_qa.ipynb',\
                            'test/pm_randoms_qa_{}.ipynb'.format(field),\
                            parameters=dict(field=field)
                            )
