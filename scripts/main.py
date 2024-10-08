import symengine
from sympy import log

from bender_wu import seriescoeff

if __name__ == "__main__":
    N = 140
    x = symengine.Symbol("x")
    fg = x + log(x - 1)
    r0 = 1.5
    prec = 100
    sf = seriescoeff.SeriesCoeff("x", r0, fg, N, prec)
    sf.taylor()
    # sf.p_coeff()
    # sf.c_coeff()
    print(sf.f_n)
