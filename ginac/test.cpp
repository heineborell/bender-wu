#include <ginac/pseries.h>
#include <iostream>
#ifdef IN_GINAC
#include "ginac.h"
#else
#include <ginac/ginac.h>
#endif
using namespace GiNaC;

void foo() {
  numeric three{3.0};
  numeric r{2, 3};
  std::cout << Pi.evalf() << '\n';
  std::cout << r.evalf() << '\n';
  std::cout << sin(r) << '\n';
}
int main() {
  symbol x{"x"};
  ex e1{sin(x)};
  ex subin{series_to_poly(e1.series(x == 0, 10))};
  Digits = 36;
  foo();
  std::cout << e1.subs(x == 7).evalf() << '\n';
  std::cout << subin.subs(x == numeric(1, 2)).evalf();

  return 0;
}
