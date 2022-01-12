# https://docs.pytest.org/en/6.2.x/
# content of test_sample.py
def inc(x):
    return x + 1


def test_answer():
    assert inc(4) == 5
