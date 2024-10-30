"""Module providing for computation of A[(l,k)]"""

import symengine

if __name__ == "__main__":

    nu = 3
    eLL = 10
    beta = symengine.symbols("beta")
    a_array = symengine.symarray(
        "a", (eLL + 1, 1 + nu + 3 * eLL)
    )  # symbolic values for A[(l,k)]
    v_array = symengine.symarray(
        "v", eLL + 3
    )  # symbolic values for potential expansion V_n
    # h_array = symengine.symarray("h", eLL + 3)  # symbolic values for function H, H_n
    h_array = [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
    e_array = symengine.symarray(
        "e", eLL + 3
    )  # symbolic values for energy levels epsilon(l+2)

    print(e_array)

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
            A[(L, k)] = 1 / (2 * (k - nu)) * A[(L, k)]

    # Computing E (energy) using the A[(l,k)] we found above
    for L in range(1, eLL + 1):
        if e_array[L + 2] != 0:
            e_array[L + 2] = -(nu + 1) * (nu + 2) * A[(L, nu + 2)]
            for n in range(1, L + 1):
                e_array[L + 2] = e_array[L + 2] + v_array[n + 2] * A[(L - n, nu - n - 2)]

        e_array[L + 2] = 1 / (2 * A[(0, nu)] * h_array[0]) * e_array[L + 2]

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

    # arrayfull = []
    # for i in range(0, eLL + 1):
    #     for j in range(0, 1 + nu + 3 * eLL):
    #         arrayfull.append((i, j))
    # a_subs_dict = dict(zip(a_array.flatten(), {k: A[k] for k in arrayfull}.values()))

    efull = symengine.symarray("e", eLL + 3)
    for i, j in enumerate(efull):
        if i % 2 != 0:
            efull[i] = 0
    efull = efull[efull != 0]

# e_subs_dict = dict(zip(efull, e_array[e_array != 0]))

# def repeat(n, x, y):
#     for _ in range(n):
#         x = x.subs(y)
#     return x

# e_array = e_array[e_array != 0]
# e_array = [
#     repeat(eLL, i, a_subs_dict)
#     for i in [repeat(eLL, j, e_subs_dict) for j in e_array]
# ]
