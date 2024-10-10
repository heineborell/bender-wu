import symengine
from symengine.lib.symengine_wrapper import solve
from sympy import log

from bender_wu import potcoeff, seriescoeff


def fg(x):
    return x + log(x - 1)


def V(r, s, L):
    return (1 - 1 / r) * ((L * (L + 1)) / r**2 + (1 - s**2) / r**3)


if __name__ == "__main__":
    N = 10
    x = symengine.Symbol("x")
    r = symengine.Symbol("r")
    # r0 = max(solve(symengine.diff(V(r, 2, 2), r), r))
    prec = 3
    R0 = solve(symengine.diff(V(x, 2, 2), x), x)
    r0 = max([sol.n(prec) for sol in iter(R0.args)])
    sf = seriescoeff.SeriesCoeff("x", 0, fg(r0 + x) - fg(r0), N, prec)
    sf.taylor()
    sf.p_coeff()
    sf.c_coeff()
    # print(sf.c_n)
    v_coef = potcoeff.PotCoeff("x", V(x, 2, 2), r0, sf.c_n, N, prec)
    v_coef.vexpand()
    # print(sf.f_n_x0[0])
    # print(vc.cn[0])
    print(v_coef.vdel)
    v_coef.vcoeff()
    print(v_coef.vc)
