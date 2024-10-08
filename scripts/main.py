import symengine
from sympy import log

from bender_wu import series_coeff

if __name__ == "__main__":
    N = 40
    x = symengine.Symbol("x")
    r0 = 1.5
    ts = series_coeff.taylor(x, x + log(x - 1), r0, N, 200, fac=False)
    tp = series_coeff.p_coeff(N, ts)
    tc = series_coeff.c_coeff(ts, N, tp)
    print(tc)
