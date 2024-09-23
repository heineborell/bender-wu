import math

# import numpy as np
# from symengine import Symbol, diff
from sympy import diff, ln, solve, symbols

# from sympy.abc import a, b, c
# from sympy.polys import QQ, RR, ring
# from sympy.polys.ring_series import rs_series

# Create the function f(x)
x, s, r0, d, L = symbols("x s r0 d L")
# fx = x + ln(s*x -1)
# fx = Function('f')(x)
x0 = 0  # expand about the point x_0
N = 100  # get first N coefficients


def V(x, s, L):
    return (1 - 1 / x) * ((L * (L + 1)) / x**2 + (1 - s**2) / x**3)


R0 = max(solve(diff(V(x, 2, 2), x), x)).n()
print(R0)


def rs(x):
    return x + ln(x - 1)


fx = rs(r0 + x) - rs(r0)


f_n = {}  # nth derivatives of f(x)
f_n_x0 = {}  # nth derivatives, but evaluated at x = x_0
f_n[1] = diff(fx, x, 1)  # differentiate f(x)
f_n_x0[1] = f_n[1].subs(
    {x: x0, r0: max(solve(diff(V(x, 2, 2), x), x)).n()}
)  # do the substitution

for i in range(2, N + 1):
    f_n[i] = diff(f_n[i - 1], x, 1)
    # it's more efficient to differentiate the previous derivative
    # once than to directly ask for the nth derivative
    f_n_x0[i] = f_n[i].subs({x: x0, r0: max(solve(diff(V(x, 2, 2), x), x)).n()}).evalf()

print(f_n_x0)


def P(N: int, Y: list[float]) -> dict[int, float]:
    # Y should contain N symbolic derivatives, where
    # Y[1] is the first derivative of f(x),
    # Y[2] is the second derivative of f(x), etc
    #
    # Returns a hash table P: (i,j) -> symbolic equation
    # where P[(i,j)] = P(i,j) from the paper
    P = {}
    for j in range(1, N + 1):
        P[(j, j)] = Y[1] ** j
        for k in range(j + 1, N + 1):
            P[(j, k)] = 0
            for m in reversed(range(1, k - j + 1)):
                P[(j, k)] = (
                    P[(j, k)]
                    + (m * j - k + j + m)
                    * Y[m + 1]
                    / math.factorial(m + 1)
                    * P[(j, k - m)]
                )
            P[(j, k)] = P[(j, k)] * 1 / (k - j) * 1 / Y[1]
    return P


# P = P(N, f_n_x0)
#
#
# b_n = {}  # Vector of pre-computed dummy variable values
# b_n[1] = 1 / f_n_x0[1]
#
# c_n = {}  # vector of Taylor series coefficients
# c_n[1] = b_n[1] / factorial(1)
#
# for n in range(2, N + 1):
#
#    b_n[n] = 0
#    for j in range(1, n):
#        b_n[n] = b_n[n] + b_n[j] / factorial(j) * P[(j, n)]
#    b_n[n] = b_n[n] * factorial(n) * -1 * b_n[1] ** n
#    c_n[n] = b_n[n] / factorial(n)
#
#
# def delta(x):
#    pdel = 0
#    for i in range(1, N + 1):
#        pdel += c_n[i] * x**i
#
#    return Poly(pdel)
#
#
# Vx = V(delta(x), 2, 2)
# V_n = {}  # nth derivatives of f(x)
# V_n_x0 = {}  # nth derivatives, but evaluated at x = x_0
# V_n[1] = diff(Vx, x, 1)  # differentiate f(x)
# V_n_x0[1] = V_n[1].subs({x: R0}).evalf()  # do the substitution
#
# for i in range(2, N + 1):
#    V_n[i] = diff(V_n[i - 1], x, 1)
#    # it's more efficient to differentiate the previous derivative
#    # once than to directly ask for the nth derivative
#    V_n_x0[i] = V_n[i].subs({x: R0}).evalf()
#
#
# print(V_n_x0)
