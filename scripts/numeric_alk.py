"""This is the numerical A[(l,k)] computation"""

import math

import symengine
from symengine.lib.symengine_wrapper import solve
from sympy import log

from bender_wu import potcoeff, seriescoeff


def fg(x):
    return x + log(x - 1)


def V(r, s, L):
    return (1 - 1 / r) * ((L * (L + 1)) / r**2 + (1 - s**2) / r**3)


if __name__ == "__main__":

    nu = 2
    eLL = 20
    N = nu + (eLL + 2)
    x = symengine.Symbol("x")
    r = symengine.Symbol("r")
    s = 2
    L = 2
    prec = 200
    R0 = solve(symengine.diff(V(x, s, L), x), x)
    r0 = max([sol.n(prec) for sol in iter(R0.args)])
    sf = seriescoeff.SeriesCoeff("x", 0, fg(r0 + x) - fg(r0), N, prec)
    sf.taylor()
    sf.p_coeff()
    sf.c_coeff()
    v_coef = potcoeff.PotCoeff("x", -V(x, 2, 2), r0, sf.c_n, N, prec)
    v_coef.vexpand()
    v_coef.vcoeff()

    h_array = dict.fromkeys(range(0, eLL + 3), 0)
    h_array[0] = 1
    v_array = v_coef.vc
    e_array = dict.fromkeys(range(0, eLL + 3), 0)
    e_array[0] = 1 / 2 * (v_array[0] / h_array[0])
    beta = -math.sqrt((h_array[0] * v_array[2] - h_array[2] * v_array[0]) / h_array[0])
    e_array[2] = -(2 * beta / h_array[0]) * (nu + 1 / 2) * 1 / 2

    A = dict.fromkeys(
        [
            (L, k)
            for L in range(0, eLL + 1)
            for k in range(-(nu + 3 * (eLL + 1)), nu + 3 * (eLL + 1))
        ],
        0,
    )

    A[(0, nu)] = 1

    # Computing A[(l,k)] from k = nu+3l down to k = nu
    for L in range(1, eLL + 1):
        for k in reversed(range(nu + 1, nu + 3 * L + 1)):  # first k from nu+3L to nu
            A[(L, k)] = -(k + 1) * (k + 2) * A[(L, k + 2)]
            for n in range(1, L + 1):
                A[(L, k)] = A[(L, k)] + v_array[n + 2] * A[(L - n, k - n - 2)]
                for m in range(0, n + 3):
                    A[(L, k)] = (
                        A[(L, k)] - h_array[m] * e_array[n + 2 - m] * A[(L - n, k - m)]
                    )

            A[(L, k)] = 1 / (2 * (k - nu) * beta) * A[(L, k)]

    # Computing E (energy) using the A[(l,k)] we found above
    for L in range(1, eLL + 1):
        e_array[L + 2] = -(nu + 1) * (nu + 2) * A[(L, nu + 2)]
        for n in range(1, L + 1):
            e_array[L + 2] = e_array[L + 2] + v_array[n + 2] * A[(L - n, nu - n - 2)]
        e_array[L + 2] = 1 / (2 * A[(0, nu)] * h_array[0]) * e_array[L + 2]

    print(A[0, 2])
    print(e_array)
