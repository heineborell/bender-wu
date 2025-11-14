#include "engine.h"
#include "helper.h"
#include <ginac/pseries.h>

int main() {
  Digits = 100;
  Timer t;
  symbol r{"ra"};
  symbol x{"x"};
  int order{30};
  ex radius{x * x * (x - 1) * (x - 1) / 2};
  std::vector<ex> result{fourierSeries(radius, x, order)};
  // std::vector<std::vector<ex>> pcoeff{pCoeff(result, order)};
  // std::vector<ex> ccoeff{cCoeff(result, pcoeff, radius, x, order)};
  std::vector<std::vector<ex>> acoeff{aCoeff(result)};
  std::vector<ex> energy{Energy(acoeff, result)};
  // ex series_expansion{radius.series(x == 0, 4)};

  // ex subin{series_to_poly(e1.series(x == 0, 10))};
  // std::cout << series_to_poly(series_expansion) << '\n';
  // printDict(acoeff);
  printArray(energy, x, false);
  // printArray(result, x, false);
  // printArray(result, x);
  std::cout << "Time elapsed: " << t.elapsed() << " seconds\n";

  return 0;
}
