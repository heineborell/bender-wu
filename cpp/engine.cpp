#include "engine.h"

std::vector<RCP<const Basic>>
fourierSeries(RCP<const Basic> func, RCP<const Symbol> var, std::size_t order) {
  std::vector<RCP<const Basic>> arr(order + 1);
  arr[0] = func;            // 0.th order
  arr[1] = func->diff(var); // 1. order
  for (std::size_t i{2}; i <= order; ++i) {
    arr[i] = arr[i - 1]->diff(var);
  }

  return arr;
}

std::vector<std::vector<RCP<const Basic>>>
inverseCoeff(std::vector<RCP<const Basic>> &arr, std::size_t order) {
  std::vector<std::vector<RCP<const Basic>>> dict(
      order + 1,
      std::vector<RCP<const Basic>>(
          order +
          1)); // in order to start the index from 1, the size is order+1
  for (std::size_t j{1}; j <= order; ++j) {
    dict[j][j] = pow(arr[1], SymEngine::integer(j));
    for (std::size_t k{j + 1}; k <= order; ++k) {
      dict[j][k] = SymEngine::integer(0);
      for (std::size_t m{k - j}; m >= 1; --m) {
        dict[j][k] = add(dict[j][k],
                         mul(mul(SymEngine::integer(m * j - k + j + m),
                                 div(arr[m + 1], SymEngine::factorial(m + 1))),
                             dict[j][k - m]));
      }
      dict[j][k] = (mul(dict[j][k], mul(SymEngine::integer(k - j), arr[1])));
    }
  }
  return dict;
}

RCP<const Basic> getPotential(RCP<const Symbol> r, const int L, const int s) {
  return mul(
      div(sub(r, integer(1)), r),
      add(div(mul(integer(L), add(integer(L), integer(1))), pow(r, integer(2))),
          div(sub(integer(1), pow(integer(s), integer(2))),
              pow(r, integer(3)))));
}
