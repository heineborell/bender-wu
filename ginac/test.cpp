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

numeric baz() {
  numeric ratio(3, 7);
  return ratio;
}
int main() {
  int t1{15};
  int t2{7};
  symbol x{"x"};
  ex e1{sin(x)};
  ex subin{series_to_poly(e1.series(x == 0, 10))};
  numeric t3{3 / numeric(2)};
  Digits = 30;
  // foo();
  // std::cout << e1.subs(x == 7) << '\n';
  // std::cout << subin.subs(x == numeric(1, 2)).evalf();
  // std::cout << t1 / t2 << '\n';
  // std::cout << baz() << '\n';
  std::cout << is_a<numeric>(t3) << '\n';

  return 0;
}
