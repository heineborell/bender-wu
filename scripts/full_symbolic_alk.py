"""Module providing for computation of A[(l,k)]"""

import math

import symengine
from symengine.lib.symengine_wrapper import solve
from sympy import log

from bender_wu import potcoeff, seriescoeff


def fg(x):
    return x + log(x - 1)


def V(r, s, ELL):
    return (1 - 1 / r) * ((ELL * (ELL + 1)) / r**2 + (1 - s**2) / r**3)


if __name__ == "__main__":

    nu = 0
    eLL = 5
    N = nu + (eLL + 2)
    x = symengine.Symbol("x")
    r = symengine.Symbol("r")
    s = 2
    ELL = 2
    prec = 10
    R0 = solve(symengine.diff(V(x, s, ELL), x), x)
    r0 = max([sol.n(prec) for sol in iter(R0.args)])
    sf = seriescoeff.SeriesCoeff("x", 0, fg(r0 + x) - fg(r0), N, prec)
    sf.taylor()
    sf.p_coeff()
    sf.c_coeff()
    v_coef = potcoeff.PotCoeff("x", -V(x, 2, 2), r0, sf.c_n, N, prec)
    v_coef.vexpand()
    v_coef.vcoeff()

    a_array = symengine.symarray(
        "a", (eLL + 1, 1 + nu + 3 * eLL)
    )  # symbolic values for A[(l,k)]
    e_array = symengine.symarray(
        "e", eLL + 3
    )  # symbolic values for energy levels epsilon(l+2)

    h_array = dict.fromkeys(range(0, eLL * 3), 0)
    h_array[0] = -1
    v_array = v_coef.vc
    e_array[0] = v_array[0] / h_array[0]
    beta = -math.sqrt((h_array[0] * v_array[2] - h_array[2] * v_array[0]) / h_array[0])
    e_array[2] = -(2 * beta / h_array[0]) * (nu + 1 / 2)

    # Odd valued energy levels are zero
    for i, j in enumerate(e_array):
        if i % 2 != 0:
            e_array[i] = 0
    print(e_array)
    print(beta)

    # First define A as a dict then populate the A[(l,k)] and epsilon(l+2)
    # values before using the master formula
    A = {}
    A[(0, nu)] = 1

    # Setting A[(l,k)] to symbolic a_l_k
    for L in range(0, eLL + 1):
        for k in range(0, 1 + nu + 3 * eLL):
            A[(L, k)] = a_array[L, k]

    # A[(0,k)] = 1 for k = nu (for normalization) and A[(0,k)]= 0 for k != nu
    for L in range(0, eLL + 1):
        for k in range(0, nu + 3 * L + 1):
            if k == nu:
                A[(0, k)] = 1
            else:
                A[(0, k)] = 0

    # A[(l,k)] = 0 for k<0
    for k in range(1, nu + 3 * (eLL + 1)):
        for L in range(0, eLL + 1):
            A[(L, -k)] = 0

    # A[(l,k)] = 0 for k> nu+3l
    for L in range(0, eLL + 1):
        for k in range(1, nu + 3 * (eLL + 1)):
            if k > nu + 3 * L:
                A[(L, k)] = 0

    # A[(L,nu)]= 0  for l >= 1 for normalization
    for L in range(1, eLL + 1):
        A[(L, nu)] = 0

    # Computing A[(l,k)] from k = nu+3l down to k = nu
    for L in range(1, eLL + 1):
        for k in reversed(range(nu + 1, nu + 3 * L + 1)):  # first k from nu+3L to nu
            for n in range(1, L + 1):
                A[(L, k)] = -(k + 1) * (k + 2) * A[(L, k + 2)]
                for n in range(1, L + 1):
                    A[(L, k)] = A[(L, k)] + v_array[n + 2] * A[(L - n, k - n - 2)]
                    for m in range(0, n + 3):
                        A[(L, k)] = (
                            A[(L, k)]
                            - h_array[m] * e_array[n + 2 - m] * A[(L - n, k - m)]
                        )
            A[(L, k)] = 1 / (2 * (k - nu) * beta) * A[(L, k)]

    # Computing E (energy) using the A[(l,k)] we found above
    for L in range(1, eLL + 1):
        sum_e1 = 0
        sum_e2 = 0
        sum_e3 = 0
        if e_array[L + 2] != 0:
            e_array[L + 2] = -(nu + 1) * (nu + 2) * A[(L, nu + 2)]

            for n in range(1, L + 1):
                sum_e1 += v_array[n + 2] * A[(L - n, nu - n - 2)]

            for n in range(1, L):
                for m in range(0, n + 3):
                    sum_e2 += h_array[m] * e_array[n + 2 - m] * A[(L - n, nu - m)]

            for n in range(1, L + 3):
                sum_e3 += h_array[n] * e_array[L + 2 - n] * A[(0, nu - n)]

        e_array[L + 2] = e_array[L + 2] + sum_e1 - sum_e2 - sum_e3
        e_array[L + 2] = 1 / (A[(0, nu)] * h_array[0]) * e_array[L + 2]

    # Computing the rest of A[(l,k)] from k = nu down to k = 0
    for L in range(1, eLL + 1):
        for k in reversed(range(0, nu)):
            for n in range(1, L + 1):
                A[(L, k)] = -(k + 1) * (k + 2) * A[(L, k + 2)]
                for n in range(1, L + 1):
                    A[(L, k)] = A[(L, k)] + v_array[n + 2] * A[(L - n, k - n - 2)]
                    for m in range(0, n + 3):
                        A[(L, k)] = (
                            A[(L, k)]
                            - h_array[m] * e_array[n + 2 - m] * A[(L - n, k - m)]
                        )
            A[(L, k)] = 1 / (2 * (k - nu) * beta) * A[(L, k)]

    # print(e_array[2])
    # print(sympify(expand(A[4, 2])))
    print(e_array)
