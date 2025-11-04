#include "engine.h"

std::vector<ex> fourierSeries(const ex &func, const symbol &var,
                              std::size_t order) {
  std::vector<ex> f_n(order + 1);
  std::vector<ex> f_n_x0(order + 1);
  f_n[0] = func; // 0.th order
  f_n_x0[0] = f_n[0].subs(var == 0);
  f_n[1] = func.diff(var); // 1. order
  f_n_x0[1] = f_n[1].subs(var == 0);
  for (std::size_t i{2}; i <= order; ++i) {
    f_n[i] = f_n[i - 1].diff(var);
    f_n_x0[i] = f_n[i].subs(var == 0);
  }

  return f_n_x0;
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
      int jj = static_cast<int>(j);
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
