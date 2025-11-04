import symengine
from symengine.lib.symengine_wrapper import solve
from sympy import log, sqrt

from bender_wu import potcoeff, seriescoeff
import time


def fg(x):
    return x + log(x + 1)


def V(r, s, L):
    return (1 - 1 / r) * ((L * (L + 1)) / r**2 + (1 - s**2) / r**3)


if __name__ == "__main__":
    start_time = time.perf_counter()
    N = 5
    x = symengine.Symbol("x")
    r = symengine.Symbol("r")
    s = 2
    L = 2
    prec = 10
    # R0 = solve(symengine.diff(V(x, s, L), x), x)
    # r0 = max([sol.n(prec) for sol in iter(R0.args)])
    sf = seriescoeff.SeriesCoeff("x", 0, fg(x), N, prec)
    sf.taylor()
    sf.p_coeff()
    sf.c_coeff()
    # print(sf.f_n_x0)
    # print(sf.P)
    # print(sf.c_n)

    end_time = time.perf_counter()
    elapsed_time = end_time - start_time
    print(f"Elapsed time is {elapsed_time:.6f} seconds")

# nu = 0
# beta = float(sqrt(v_coef.vc[2]))
# epsilon = {}
# epsilon[2] = 1 / 2 * beta
# eLL = 8
# A = {}
# A[(0, nu)] = 1
# A[(0, 0)] = 1
# for L in range(0, eLL):
#     for k in range(0, nu + 3 * L + 1):
#         if k > nu:
#             A[(0, k)] = 0
# for k in range(1, nu + 3 * eLL):
#     for L in range(0, nu + 3 * eLL):
#         A[(L, -k)] = 0
#
# for L in range(0, eLL):
#     for k in range(1, nu + 3 * eLL):
#         if k > nu + 3 * L:
#             A[(L, k)] = 0
#
#
# for L in range(1, eLL):
#     A[(L, nu)] = 0
#     for k in reversed(range(nu + 1, nu + 3 * L + 1)):
#         for n in range(1, L + 1):
#             A[(L, k)] = -(k + 1) * (k + 2) * A[(L, k + 2)]
#             for n in range(1, L + 1):
#                 A[(L, k)] = A[(L, k)] + float(v_coef.vc[n + 2]) * A[(L - n, k - n - 2)]
#
#         A[(L, k)] = 1 / (2 * (k - nu) * beta) * A[(L, k)]
#
#
# for L in range(1, eLL):
#     epsilon[L + 2] = -(nu + 1) * (nu + 2) * A[(L, nu + 2)]
#     for n in range(1, L + 1):
#         epsilon[L + 2] = epsilon[L + 2] + float(v_coef.vc[n + 2]) * A[(L - n, nu - n - 2)]
#
#     epsilon[L + 2] = 1 / A[(0, nu)] * epsilon[L + 2]
# print(epsilon)
