import pandas as pd


def test_func():
    # test passes on pandas.__version__ == u'0.18.0'
    df = pd.DataFrame([[1, 2, 1, 1],
                       [2, 2, 2, 2],
                       [1, 2, 3, 1],
                       [3, 2, 4, 2]], columns=['A', 'B', 'C', 'D'])

    assert all(df.duplicated('A', keep=False) == [True, False, True, False])
    assert all(df.duplicated('B', keep=False) == [True, True, True, True])
    assert all(df.duplicated('C', keep=False) == [False, False, False, False])
    assert all(df.duplicated('D', keep=False) == [True, True, True, True])

    assert all(df.duplicated('A') == [False, False, True, False])
    assert all(df.duplicated('B') == [False, True, True, True])
    assert all(df.duplicated('C') == [False, False, False, False])
    assert all(df.duplicated('D') == [False, False, True, True])
