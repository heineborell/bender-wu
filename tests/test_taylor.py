from sympy import symbols

from bender_wu.taylor import taylor


def test_taylor():
    x = symbols("x")
    want = {1: 3, 2: 6, 3: 6, 4: 0}
    got = taylor(x, x**3, 1, 4)
    assert got == want
