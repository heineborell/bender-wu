from sympy import diff, symbols


def taylor(var, func, val: float, N: int):
    var = symbols(str(var))
    f_n = {}  # nth derivatives of f(x)
    f_n_x0 = {}  # nth derivatives, but evaluated at x = x_0
    f_n[1] = diff(func, var, 1)  # differentiate f(x)
    f_n_x0[1] = f_n[1].subs(var, val)  # do the substitution

    for i in range(2, N + 1):
        f_n[i] = diff(f_n[i - 1], var, 1)
        # it's more efficient to differentiate the previous derivative
        # once than to directly ask for the nth derivative
        f_n_x0[i] = f_n[i].subs(var, val)

    return f_n_x0
