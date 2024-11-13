import symengine
from sympy.polys.domains import RR
from sympy.polys.ring_series import rs_pow
from sympy.polys.rings import ring

from bender_wu import seriescoeff


class PotCoeff:
    def __init__(self, var_pot, pot, val_pot, cn, N, prec) -> None:
        self.var_pot = var_pot
        self.pot = pot
        self.val_pot = val_pot
        self.N = N
        self.prec = prec
        self.cn = cn
        self.pdel = 0
        self.vdel = 0

        sc = seriescoeff.SeriesCoeff(
            self.var_pot, self.val_pot, self.pot, self.N, self.prec
        )
        sc.taylor(fac=True)
        self.fn = sc.f_n_x0

        R, self.z = ring("z", RR)
        for i in range(0, self.N + 1):
            self.pdel += self.cn[i] * self.z**i

    def vexpand(self):
        for i in range(0, self.N + 1):
            self.vdel += self.fn[i] * rs_pow(self.pdel, i, self.z, self.N)

    def vcoeff(self):
        self.vc = {}
        for i in range(0, self.N + 1):
            self.vc[i] = symengine.Float(self.vdel.coeff(self.z**i))
