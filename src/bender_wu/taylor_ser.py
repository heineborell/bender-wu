import math

import symengine


def taylor(var, func, val: float, N: int, repl=True) -> dict[int, float]:

    f_n = {}  # nth derivatives of f(x)
    f_n_x0 = {}  # nth derivatives, but evaluated at x = x_0
    f_n[1] = symengine.diff(func, var, 1)  # differentiate f(x)
    f_n_x0[1] = f_n[1].subs(var, val).evalf()  # do the substitution

    for i in range(2, N + 1):
        f_n[i] = symengine.diff(f_n[i - 1], var, 1)
        # it's more efficient to differentiate the previous derivative
        # once than to directly ask for the nth derivative
        f_n_x0[i] = f_n[i].subs(var, val).evalf()

    if repl is True:
        return f_n_x0
    else:
        return f_n


def p_coeff(N: int, Y: dict) -> dict[int, float]:
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


def c_coeff(f_n_x0, N, P):
    b_n = {}  # Vector of pre-computed dummy variable values

    b_n[1] = 1 / f_n_x0[1]

    c_n = {}  # vector of Taylor series coefficients
    c_n[1] = b_n[1] / math.factorial(1)

    for n in range(2, N + 1):

        b_n[n] = 0
        for j in range(1, n):
            b_n[n] = b_n[n] + b_n[j] / math.factorial(j) * P[(j, n)]
            b_n[n] = b_n[n] * math.factorial(n) * -1 * b_n[1] ** n
            c_n[n] = b_n[n] / math.factorial(n)
    return c_n
