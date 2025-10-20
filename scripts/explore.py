import symengine
from symengine.lib.symengine_wrapper import solve
from sympy import log
import time
import math

from bender_wu import potcoeff, seriescoeff


def fg(x):
    return x + log(x - 1)


def V(r, s, L):
    return (1 - 1 / r) * ((L * (L + 1)) / r**2 + (1 - s**2) / r**3)


if __name__ == "__main__":
    start_time = time.perf_counter()
    N = 3
    x = symengine.Symbol("x")
    r = symengine.Symbol("r")
    s = 2
    L = 2
    prec = 10
    R0 = solve(symengine.diff(V(x, s, L), x), x)
    r0 = max([sol.n(prec) for sol in iter(R0.args)])
    print(R0)
    print(r0)
    sf = seriescoeff.SeriesCoeff("x", 0, fg(r0 + x) - fg(r0), N, prec)
    sf.taylor(fac=True)
    sf.p_coeff()
    # sf.c_coeff()
    print(sf.f_n)
    # print((sf.P))
    # print(len(sf.P))

    end_time = time.perf_counter()
    elapsed_time = end_time - start_time
    print(f"Elapsed time is {elapsed_time:.6f} seconds")
