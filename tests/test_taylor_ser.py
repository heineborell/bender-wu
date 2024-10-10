from sympy import symbols

from bender_wu import seriescoeff


def test_taylor():
    x = symbols("x")
    want = {0: 1.0, 1: 3.0, 2: 6.0, 3: 6.0, 4: 0.0}
    ts = seriescoeff.SeriesCoeff("x", 1, x**3, 4, 2)
    ts.taylor()
    got = ts.f_n_x0
    assert got == want


def test_taylor_fac():
    x = symbols("x")
    want = {0: 1.0, 1: 3.0 / 1, 2: 6.0 / 2, 3: 6.0 / 6, 4: 0.0 / 12}
    ts = seriescoeff.SeriesCoeff("x", 1, x**3, 4, 2)
    ts.taylor(fac=True)
    got = ts.f_n_x0
    assert got == want
