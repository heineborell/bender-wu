import symengine

from bender_wu import taylor_ser

if __name__ == "__main__":
    x = symengine.Symbol("x")
    tc = taylor_ser.taylor(x, x**4, 1, 5)
    print(taylor_ser.p_coeff(3, tc))
    print(tc)
