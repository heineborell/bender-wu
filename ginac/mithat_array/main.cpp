#include "engine.h"
#include "helper.h"
#include <ginac/pseries.h>

int main() {
  Digits = 20;
  Timer t;
  symbol r{"ra"};
  symbol x{"x"};
  ex potential{(x * x) + (x * x * x * x)};
  ex omega{sqrt(potential.diff(x, 2)).subs(x == 0)};
  std::array<ex, 2 * expansionOrder + 1> result{vSeries(potential, x)};
  std::vector<ex> ecoeff{energy(result, omega)};

  printArray(ecoeff, true);
  std::cout << "Time elapsed: " << t.elapsed() << " seconds\n";

  return 0;
}
