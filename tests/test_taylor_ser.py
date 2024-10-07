from sympy import symbols

from bender_wu import taylor_ser


def test_taylor():
    x = symbols("x")
    want = {1: 3.0, 2: 6.0, 3: 6.0, 4: 0.0}
    got = taylor_ser.taylor(x, x**3, 1, 4, fac=False)
    assert got == want


def test_taylor_fac():
    x = symbols("x")
    want = {1: 3.0 / 1, 2: 6.0 / 2, 3: 6.0 / 6, 4: 0.0 / 12}
    got = taylor_ser.taylor(x, x**3, 1, 4, fac=True)
    assert got == want
