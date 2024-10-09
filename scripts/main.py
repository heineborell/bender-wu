import symengine
from symengine.lib.symengine_wrapper import solve
from sympy import log

from bender_wu import seriescoeff


def fg(x):
    return x + log(x - 1)


def V(r, s, L):
    return (1 - 1 / r) * ((L * (L + 1)) / r**2 + (1 - s**2) / r**3)


if __name__ == "__main__":
    N = 300
    x = symengine.Symbol("x")
    r = symengine.Symbol("r")
    # r0 = max(solve(symengine.diff(V(r, 2, 2), r), r))
    prec = 200
    R0 = solve(symengine.diff(V(x, 2, 2), x), x)
    r0 = max([sol.n(prec) for sol in iter(R0.args)])
    sf = seriescoeff.SeriesCoeff("x", 0, fg(r0 + x) - fg(r0), N, prec)
    sf.taylor()
    sf.p_coeff()
    sf.c_coeff()
    print(sf.c_n[100])
