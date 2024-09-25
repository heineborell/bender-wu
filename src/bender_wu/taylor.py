from sympy import diff, symbols


def taylor(func, val, N):

    f_n = {}  # nth derivatives of f(x)
    f_n_x0 = {}  # nth derivatives, but evaluated at x = x_0
    f_n[1] = diff(func, x, 1)  # differentiate f(x)
    f_n_x0[1] = f_n[1].subs(x, val)  # do the substitution

    for i in range(2, N + 1):
        f_n[i] = diff(f_n[i - 1], x, 1)
        # it's more efficient to differentiate the previous derivative
        # once than to directly ask for the nth derivative
        f_n_x0[i] = f_n[i].subs(x, val).evalf()

    return f_n_x0


if __name__ == "__main__":

    x = symbols("x")

    print(taylor(x**3, 1, 4))
