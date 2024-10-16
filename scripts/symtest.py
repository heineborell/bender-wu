import symengine

if __name__ == "__main__":

    A = {}
    nu = 0
    eLL = 5
    beta = symengine.symbols("beta")
    aarray = symengine.symarray("a", (eLL + 1, 1 + nu + 3 * eLL))
    varray = symengine.symarray("v", 20)
    # harray = symengine.symarray("h", 20)
    harray = [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
    earray = symengine.symarray("e", 20)

    for i, j in enumerate(earray):
        if i % 2 != 0:
            earray[i] = 0

    for L in range(0, eLL + 1):
        for k in range(0, nu + 3 * eLL):
            A[(L, k)] = aarray[L, k]
    for L in range(0, eLL + 1):
        for k in range(0, nu + 3 * L + 1):
            if k > nu:
                A[(0, k)] = 0
    for k in range(1, nu + 3 * eLL):
        for L in range(0, eLL + 1):
            A[(L, -k)] = 0
    for L in range(0, eLL + 1):
        for k in range(1, nu + 3 * (eLL + 1)):
            if k > nu + 3 * L:
                A[(L, k)] = 0

    A[(0, nu)] = 1
    for L in range(1, eLL + 1):
        A[(L, nu)] = 0
        for k in reversed(range(nu + 1, nu + 3 * L + 1)):
            for n in range(1, L + 1):
                A[(L, k)] = -(k + 1) * (k + 2) * A[(L, k + 2)]
                for n in range(1, L + 1):
                    A[(L, k)] = A[(L, k)] + varray[n + 2] * A[(L - n, k - n - 2)]
                    for m in range(0, n + 3):
                        A[(L, k)] = (
                            A[(L, k)] - harray[m] * earray[n + 2 - m] * A[(L - n, k - m)]
                        )
            A[(L, k)] = 1 / (2 * (k - nu) * beta) * A[(L, k)]

    for L in range(1, eLL + 1):
        if earray[L + 2] != 0:
            earray[L + 2] = -(nu + 1) * (nu + 2) * A[(L, nu + 2)]
            for n in range(1, L + 1):
                earray[L + 2] = earray[L + 2] + varray[n + 2] * A[(L - n, nu - n - 2)]

        earray[L + 2] = 1 / (2 * A[(0, nu)] * harray[0]) * earray[L + 2]

    print(earray)
