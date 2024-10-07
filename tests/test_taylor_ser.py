from sympy import symbols

from bender_wu.taylor_ser import taylor


def test_taylor():
    x = symbols("x")
    want = {1: 3.0, 2: 6.0, 3: 6.0, 4: 0.0}
    got = taylor(x, x**3, 1, 4)
    assert got == want
