#include "engine.h"
#include "helper.h"
#include <ginac/pseries.h>

int main() {
  Timer t;
  symbol r{"ra"};
  symbol x{"x"};
  ex radius{log(1 + x)};
  std::vector<ex> result{fourierSeries(radius, x, 3)};
  std::vector<std::vector<ex>> coeff{inverseCoeff(result, 3)};
  // ex series_expansion{radius.series(x == 0, 4)};

  // ex subin{series_to_poly(e1.series(x == 0, 10))};
  // Digits = 36;
  // printArray(result);
  // std::cout << series_to_poly(series_expansion) << '\n';
  printDict(coeff);
  std::cout << "Time elapsed: " << t.elapsed() << " seconds\n";

  return 0;
}
