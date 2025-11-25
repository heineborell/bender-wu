#include "engine.h"
#include <algorithm>
#include <cstddef>

std::vector<ex> fourierSeries(const ex &func, const symbol &var, int order) {
  std::vector<ex> f_n(static_cast<std::size_t>(order + 1));
  std::vector<ex> f_n_x0(static_cast<std::size_t>(order + 1));

  f_n[0] = func; // 0.th order
  f_n_x0[0] = f_n[0].subs(var == 0);
  f_n[1] = func.diff(var); // 1. order
  f_n_x0[1] = f_n[1].subs(var == 0);

  for (int i{2}; i <= order; ++i) {
    f_n.data()[i] = f_n.data()[i - 1].diff(var);
    f_n_x0.data()[i] = f_n.data()[i].subs(var == 0);
  }

  return f_n_x0;
}

std::vector<ex> vSeries(const ex &func, const symbol &var, int order) {
  std::vector<ex> f_n(static_cast<std::size_t>(order + 1));
  std::vector<ex> f_n_x0(static_cast<std::size_t>(order + 1));

  f_n[0] = func; // vn begins at 1 order with third derivative i.e.
                 // D[f,{x,i+2}]/(i+2)! {i,1,lmax}
  f_n_x0[0] = f_n[0].subs(var == 0);
  f_n[1] = func.diff(var, 3) / factorial(3); // 3. order
  f_n_x0[1] = f_n[1].subs(var == 0);

  for (int i{2}; i <= order; ++i) {
    f_n.data()[i] = f_n.data()[i - 1].diff(var) / (i + 2);
    f_n_x0.data()[i] = f_n.data()[i].subs(var == 0);
  }

  return f_n_x0;
}

void evenCheck(const std::vector<ex> &f_n_x0) {
  ex sum{};
  for (int i{1}; i <= lmax; i += 2) {
    sum += f_n_x0.data()[i] * f_n_x0.data()[i];
  }
  if (sum == 0) {
    lstep = 2;
    even = 1;
    std::cout << "lstep is 2" << '\n';
  }
  if (sum != 0) {
    lstep = 1;
    even = 0;
    std::cout << "lstep is 1" << '\n';
  }
}

std::vector<std::vector<ex>> pCoeff(const std::vector<ex> &f_n_x0,
                                    std::size_t order) {
  std::vector<std::vector<ex>> dict(
      order + 1,
      std::vector<ex>(
          order +
          1)); // in order to start the index from 1, the size is order+1

  for (std::size_t j{1}; j <= order; ++j) {
    dict[j][j] = pow(f_n_x0[1], j);
    for (std::size_t k{j + 1}; k <= order; ++k) {
      dict[j][k] = 0; // start all P[j,k] with zero, k= j+1,j+2,...
                      // then populate P[j,k] starting from top

      int jj = static_cast<int>(
          j); // cast these to ints so that you don't get wrapped around results
      int kk = static_cast<int>(k);

      for (std::size_t m{k - j}; m >= 1; --m) {
        int mm = static_cast<int>(m);
        dict[j][k] = dict[j][k] + (mm * jj - kk + jj + mm) * f_n_x0[m + 1] /
                                      factorial(mm + 1) * dict[j][k - m];
      }
      dict[j][k] = dict[j][k] * 1 / (kk - jj) * 1 / f_n_x0[1];
    }
  }
  return dict;
}

std::vector<ex> cCoeff(const std::vector<ex> &f_n_x0,
                       const std::vector<std::vector<ex>> &dict, const ex &func,
                       const symbol &var, std::size_t order) {

  std::vector<ex> b_n(order + 1);
  std::vector<ex> c_n(order + 1);

  b_n[1] = 1 / f_n_x0[1];
  c_n[0] = func.subs(var == 0);
  c_n[1] = b_n[1] / factorial(1);

  for (std::size_t n{2}; n <= order; ++n) {
    b_n[n] = 0;
    for (std::size_t j{1}; j <= n - 1; ++j) {
      b_n[n] = b_n[n] + b_n[j] / factorial(j) * dict[j][n];
    }
    b_n[n] = b_n[n] * factorial(n) * (-1) * pow(b_n[1], n);
    c_n[n] = b_n[n] / factorial(n);
  }

  return c_n;
}

std::vector<ex> energy(const std::vector<ex> &f_n_x0, ex &omega) {
  evenCheck(f_n_x0);

  std::vector<std::vector<ex>> A(lmax + 1, std::vector<ex>(nu + 3 * lmax + 3));
  A[0][nu] = 1; // Normalization

  std::vector<ex> E(lmax / 2 + 1);      // Energy
  E[0] = omega * (nu + numeric(1) / 2); // 0.th level
  ex sum{0};

  for (int l{0}; l <= lmax; l += lstep) {
    std::cout << "Computing order " << l / 2 - 1 << '\n';
    if (l > 0) {
      // Compute A[l][k] for k>nu, l>0
      A.data()[l][nu] = 0;
      for (int k{nu + 3 * l}; k > nu; --k) {
        A.data()[0].data()[k] = 0; // A[0][k]=0 as Hermite is bounded with nu
        sum = numeric(k + 2) * numeric(k + 1) * A.data()[l].data()[k + 2];

        for (int n{1};
             n <=
             std::min({k - 2, l, static_cast<int>(std::ssize(f_n_x0) - 1)});
             ++n) {
          sum += -numeric(2) * f_n_x0.data()[n] *
                 A.data()[l - n].data()[k - n - 2];
        }

        for (int n{1}; n <= l / 2; ++n) {
          sum += 2 * E.data()[n] * A.data()[l - 2 * n].data()[k];
        }
        A.data()[l].data()[k] = evalf(sum / (omega * numeric(2 * (k - nu))));

        // Compute E[l/2]
        if (l % 2 == 0) {
          sum = -numeric(nu + 2) * numeric(nu + 1) *
                A.data()[l].data()[nu + 2] / 2;
          for (int n{1}; n <= std::min(l, nu - 2); ++n) {
            sum += f_n_x0.data()[n] * A.data()[l - n].data()[nu - n - 2];
          }
          E.data()[l / 2] = evalf(sum);
        }
      }
    }

    // Compute A[l][k] for k<nu-1 for any l>=0
    for (int k{nu - 1}; k >= 0; --k) {
      sum = numeric(k + 2) * numeric(k + 1) * A.data()[l].data()[k + 2];
      for (int n{1}; n <= std::min({k - 2, l}); ++n) {
        sum +=
            -numeric(2) * f_n_x0.data()[n] * A.data()[l - n].data()[k - n - 2];
      }
      for (int n{1}; n <= l / 2; ++n) {
        sum += 2 * E.data()[n] * A.data()[l - 2 * n].data()[k];
      }
      A.data()[l].data()[k] = evalf(sum / (omega * numeric(2 * (k - nu))));
    }
  }
  return E;
}
