import papermill as pm

# https://docs.pytest.org/en/6.2.x/
def inc(x):
    return x + 1

def test_answer():
    assert inc(4) == 5

def test_notebooks():
    fields = ['G9', 'G12', 'G15']
    
    for field in ['G9']:    
        pm.execute_notebook('docs/nb/randoms_qa.ipynb',\
                            'test/pm_randoms_qa.ipynb',\
                            parameters=dict(field=field)
                            )
