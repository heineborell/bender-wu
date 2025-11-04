#include "engine.h"
#include "helper.h"
#include <ginac/pseries.h>

int main() {
  Digits = 10;
  Timer t;
  symbol r{"ra"};
  symbol x{"x"};
  ex radius{x + log(1 + x)};
  std::vector<ex> result{fourierSeries(radius, x, 5)};
  std::vector<std::vector<ex>> pcoeff{pCoeff(result, 5)};
  std::vector<ex> ccoeff{cCoeff(result, pcoeff, radius, x, 5)};
  // ex series_expansion{radius.series(x == 0, 4)};

  // ex subin{series_to_poly(e1.series(x == 0, 10))};
  // Digits = 36;
  // printArray(result);
  // std::cout << series_to_poly(series_expansion) << '\n';
  // printDict(pcoeff);
  // printArray(ccoeff, x);
  // printArray(result, x);
  std::cout << "Time elapsed: " << t.elapsed() << " seconds\n";

  return 0;
}
