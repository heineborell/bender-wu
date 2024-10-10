import math

import symengine


class SeriesCoeff:
    def __init__(self, var: str, val: float, func, N: int, prec: int) -> None:
        self.var = symengine.Symbol(var)
        self.val = val
        self.func = func
        self.N = N
        self.prec = prec
        self.f_n = {}  # nth derivatives of f(x)
        self.f_n_x0 = {}  # nth derivatives, but evaluated at x = x_0
        self.f_n[0] = self.func
        self.f_n[1] = symengine.diff(self.func, self.var, 1)  # differentiate f(x)
        self.f_n_x0[0] = (
            self.f_n[0].subs(self.var, self.val).n(self.prec)
        )  # do the substitution

        self.f_n_x0[1] = (
            self.f_n[1].subs(self.var, self.val).n(self.prec)
        )  # do the substitution
        self.P = {}
        self.b_n = {}  # Vector of pre-computed dummy variable values
        self.c_n = {}

    def taylor(self, fac=False):
        for i in range(2, self.N + 1):
            self.f_n[i] = symengine.diff(self.f_n[i - 1], self.var, 1) / (
                i if fac else 1  # differentiated the previous f as it is faster
            )
            self.f_n_x0[i] = self.f_n[i].subs(self.var, self.val).n(self.prec)

    def p_coeff(self):
        # Y should contain N symbolic derivatives, where
        # Y[1] is the first derivative of f(x),
        # Y[2] is the second derivative of f(x), etc
        #
        # Returns a hash table P: (i,j) -> symbolic equation
        # where P[(i,j)] = P(i,j) from the paper
        for j in range(1, self.N + 1):
            self.P[(j, j)] = self.f_n_x0[1] ** j
            for k in range(j + 1, self.N + 1):
                self.P[(j, k)] = 0
                for m in reversed(range(1, k - j + 1)):
                    self.P[(j, k)] = (
                        self.P[(j, k)]
                        + (m * j - k + j + m)
                        * self.f_n_x0[m + 1]
                        / math.factorial(m + 1)
                        * self.P[(j, k - m)]
                    )
                self.P[(j, k)] = self.P[(j, k)] * 1 / (k - j) * 1 / self.f_n_x0[1]

    def c_coeff(self):
        self.b_n[1] = 1 / self.f_n_x0[1]
        self.c_n[0] = self.func.subs(self.var, self.val).n(self.prec)
        self.c_n[1] = self.b_n[1] / math.factorial(1)
        for n in range(2, self.N + 1):

            self.b_n[n] = 0
            for j in range(1, n):
                self.b_n[n] = (
                    self.b_n[n] + self.b_n[j] / math.factorial(j) * self.P[(j, n)]
                )

            self.b_n[n] = self.b_n[n] * math.factorial(n) * -1 * self.b_n[1] ** n
            self.c_n[n] = self.b_n[n] / math.factorial(n)
