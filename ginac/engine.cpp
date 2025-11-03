#include "engine.h"

std::vector<ex> fourierSeries(const ex &func, symbol &var, std::size_t order) {
  std::vector<ex> arr(order + 1);
  arr[0] = func;           // 0.th order
  arr[1] = func.diff(var); // 1. order
  for (std::size_t i{2}; i <= order; ++i) {
    arr[i] = arr[i - 1].diff(var);
  }

  return arr;
}

std::vector<std::vector<ex>> pCoeff(std::vector<ex> &arr, std::size_t order) {
  std::vector<std::vector<ex>> dict(
      order + 1,
      std::vector<ex>(
          order +
          1)); // in order to start the index from 1, the size is order+1
  for (std::size_t j{1}; j <= order; ++j) {
    dict[j][j] = pow(arr[1], j);
    for (std::size_t k{j + 1}; k <= order; ++k) {
      dict[j][k] = 0; // start all P[j,k] with zero, k= j+1,j+2,...
      // then populate P[j,k] starting from top
      for (std::size_t m{k - j}; m >= 1; --m) {
        dict[j][k] = dict[j][k] + (m * j - k + j + m) * arr[m + 1] /
                                      factorial(m + 1) * dict[j][k - m];
      }
      dict[j][k] = dict[j][k] * 1 / (k - j) * 1 / arr[1];
    }
  }
  return dict;
}
