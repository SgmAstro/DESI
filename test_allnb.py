import papermill as pm

# https://docs.pytest.org/en/6.2.x/
def inc(x):
    return x + 1

def test_answer():
    assert inc(4) == 5

def test_notebooks():
    fields = ['G9', 'G12', 'G15']
    
    for field in fields:    
        pm.execute_notebook('docs/nb/randoms_qa.ipynb',\
                            'test/pm_randoms_qa_{}.ipynb'.format(field),\
                            parameters=dict(field=field)
                            )
