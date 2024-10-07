import symengine
from sympy import log

from bender_wu import taylor_ser

if __name__ == "__main__":
    N = 30
    x = symengine.Symbol("x")
    r0 = 1.5
    ts = taylor_ser.taylor(x, x + log(x - 1), 1.5, N)
    tp = taylor_ser.p_coeff(N, ts)
    tc = taylor_ser.c_coeff(ts, N, tp)
    tss = symengine.series((1 - 1 / x) * ((6 / x**2) - 3 / x**3), x, n=10)
    print(tss)
