#include "engine.h"
#include "helper.h"
#include <ginac/pseries.h>

int main() {
  Digits = 10;
  Timer t;
  symbol r{"ra"};
  symbol x{"x"};
  int order{20};
  ex radius{x * x + x * x * x * x};
  ex omega{sqrt(radius.diff(x, 2).subs(x == 0))};
  ex potential{radius / (omega * omega)};
  std::vector<ex> result{fourierSeries(potential, x, order)};
  // std::vector<std::vector<ex>> pcoeff{pCoeff(result, order)};
  // std::vector<ex> ccoeff{cCoeff(result, pcoeff, radius, x, order)};
  std::vector<std::vector<ex>> acoeff{aCoeff(result, omega)};
  std::vector<ex> energy{Energy(acoeff, result, omega)};
  // ex series_expansion{radius.series(x == 0, 4)};

  // ex subin{series_to_poly(e1.series(x == 0, 10))};
  // std::cout << series_to_poly(series_expansion) << '\n';
  //
  printDict(acoeff, false);
  printArray(energy, x, true);
  // printArray(result, x, false);
  // printArray(result, x);
  std::cout << "Time elapsed: " << t.elapsed() << " seconds\n";

  return 0;
}
