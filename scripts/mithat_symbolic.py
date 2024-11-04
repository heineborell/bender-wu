"""This is the numerical A[(l,k)] computation"""

import symengine
from sympy.polys.rings import sympify

if __name__ == "__main__":

    nu = 0
    eLL = 4
    N = nu + (eLL + 2)
    prec = 10
    beta = symengine.Symbol(r"β")
    a_array = symengine.symarray(
        "a", (eLL + 1, 1 + nu + 3 * eLL)
    )  # symbolic values for A[(l,k)]
    v_array = symengine.symarray(
        "v", eLL + 3
    )  # symbolic values for potential expansion V_n
    # h_array = symengine.symarray("h", eLL + 3)  # symbolic values for function H, H_n
    e_array = symengine.symarray(
        "e", eLL + 3
    )  # symbolic values for energy levels epsilon(l+2)

    h_array = dict.fromkeys(range(0, eLL + 3), 0)
    h_array[0] = 1
    v_array[1] = 0
    v_array[3] = 0
    v_array[4] = 0
    e_array[0] = beta * (nu + 1 / 2)

    # First define A as a dict then populate the A[(l,k)] and epsilon(l+2)
    # values before using the master formula
    A = {}
    A[(0, nu)] = 1

    # Odd valued energy levels are zero
    for i, j in enumerate(e_array):
        if i % 2 != 0:
            e_array[i] = 0

    # Setting A[(l,k)] to symbolic a_l_k
    for L in range(0, eLL + 1):
        for k in range(0, 1 + nu + 3 * eLL):
            A[(L, k)] = a_array[L, k]

    # A[(0,k)] = 1 for k = nu (for normalization) and A[(0,k)]= 0 for k != nu
    for L in range(0, eLL + 1):
        for k in range(0, nu + 3 * L + 1):
            if k == nu:
                A[(0, k)] = 1
            elif k > nu:
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
    sum_a1 = dict.fromkeys(
        [
            (L, k)
            for L in range(0, eLL + 1)
            for k in range(-(nu + 3 * (eLL + 1)), nu + 3 * (eLL + 1))
        ],
        0,
    )
    for L in range(1, eLL + 1):
        for k in reversed(range(nu + 1, nu + 3 * L + 1)):  # first k from nu+3L to nu
            for n in range(1, L + 1):
                sum_a1[(L, k)] += -2 * v_array[n] * A[(L - n, k - n - 2)]

        # A[(L, k)] = A[(L, k)] + (k + 1) * (k + 2) * A[(L, k + 2)]
        # A[(L, k)] = 1 / (2 * (k - nu) * beta) * A[(L, k)]

    sum_a2 = dict.fromkeys(
        [
            (L, k)
            for L in range(0, eLL + 1)
            for k in range(-(nu + 3 * (eLL + 1)), nu + 3 * (eLL + 1))
        ],
        0,
    )
    for L in range(1, eLL + 1):
        for k in reversed(range(nu + 1, nu + 3 * L + 1)):  # first k from nu+3L to nu
            for n in range(1, L):
                sum_a2[(L, k)] += 2 * e_array[n] * A[(L - n, k)]

    for L in range(1, eLL + 1):
        for k in reversed(range(nu + 1, nu + 3 * L + 1)):  # first k from nu+3L to nu
            A[(L, k)] = sum_a1[(L, k)] + (k + 1) * (k + 2) * A[(L, k + 2)]
            A[(L, k)] = (1 / (2 * (k - nu) * beta)) * A[(L, k)]
    print(sympify(A[(4, 2)].subs({"a_2_6": 0})))

    print(sum_a1[(2, 4)])

    # for L in range(1, eLL + 1):
    #    for k in reversed(range(nu + 1, nu + 3 * L + 1)):  # first k from nu+3L to nu
    #        A[(L, k)] = 1 / (2 * (k - nu) * beta) * A[(L, k)]

    # Computing A[(l,k)] from k = nu+3l down to k = nu
    # for L in range(1, eLL + 1):
    #    for k in reversed(range(nu + 1, nu + 3 * L + 1)):  # first k from nu+3L to nu
    #        A[(L, k)] = (k + 1) * (k + 2) * A[(L, k + 2)]

    # for L in range(1, eLL + 1):
    #    for k in reversed(range(nu + 1, nu + 3 * L + 1)):  # first k from nu+3L to nu
    #        for n in range(1, L + 1):
    #            A[(L, k)] = A[(L, k)] - 2 * v_array[n] * A[(L - n, k - n - 2)]

    # Computing E (energy) using the A[(l,k)] we found above
    # for L in range(1, eLL + 1):
    #    e_array[L] = -1 / 2 * (nu + 1) * (nu + 2) * A[(L, nu + 2)]
    #    for n in range(1, L + 1):
    #        e_array[L] = e_array[L] + v_array[n] * A[(L - n, nu - n - 2)]

    # for L in range(1, eLL + 1):
    #    print("for L", L)
    #    for k in reversed(range(nu + 1, nu + 3 * L + 1)):  # first k from nu+3L to nu
    #        print("for k", k)
    #        for n in range(1, L + 1):
    #            print("for n", n)
