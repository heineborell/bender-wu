import symengine
from sympy import log

from bender_wu import taylor_ser

if __name__ == "__main__":
    N = 300
    x = symengine.Symbol("x")
    ts = taylor_ser.taylor(x, x * log(x - 1), 3, N)
    tp = taylor_ser.p_coeff(N, ts)
    print(ts)
