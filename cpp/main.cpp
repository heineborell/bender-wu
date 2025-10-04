#include <chrono>

#include <cstddef>
#include <symengine/add.h>
#include <symengine/basic.h>
#include <symengine/functions.h>
#include <symengine/integer.h>
#include <symengine/pow.h>
#include <symengine/rational.h>
#include <symengine/series.h>
#include <symengine/symengine_casts.h>
#include <symengine/symengine_rcp.h>

using SymEngine::add;
using SymEngine::Basic;
using SymEngine::integer;
using SymEngine::log;
using SymEngine::Pow;
using SymEngine::pow;
using SymEngine::RCP; // these are reference counted pointers so its working
                      // like python under the hood.
using SymEngine::Symbol;
using SymEngine::symbol;

class Timer {
private:
  // Type aliases to make accessing nested type easier
  using Clock = std::chrono::steady_clock;
  using Second = std::chrono::duration<double, std::ratio<1>>;

  std::chrono::time_point<Clock> m_beg{Clock::now()};

public:
  void reset() { m_beg = Clock::now(); }

  double elapsed() const {
    return std::chrono::duration_cast<Second>(Clock::now() - m_beg).count();
  }
};

void printDict(std::vector<std::vector<RCP<const Basic>>> &dict) {
  std::size_t rows{dict.size()};
  std::size_t columns{dict[0].size()};
  std::cout << "the dictionary is " << rows - 1 << 'x' << columns - 1 << '\n';
  for (std::size_t i{1}; i < rows; ++i) {
    for (std::size_t j{1}; j < columns; ++j) {
      std::cout << *dict[i][j];
    }
  }
}
void printArray(std::vector<RCP<const Basic>> &arr) {
  for (RCP<const Basic> item : arr) {
    std::cout << *item << '\n';
  }
}

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

// def p_coeff(self):
//     for j in range(1, self.N + 1):
//         self.P[(j, j)] = self.f_n_x0[1] ** j
//         for k in range(j + 1, self.N + 1):
//             self.P[(j, k)] = 0
//             for m in reversed(range(1, k - j + 1)):
//                 self.P[(j, k)] = (
//                     self.P[(j, k)]
//                     + (m * j - k + j + m)
//                     * self.f_n_x0[m + 1]
//                     / math.factorial(m + 1)
//                     * self.P[(j, k - m)]
//                 )
//             self.P[(j, k)] = self.P[(j, k)] * 1 / (k - j) * 1 /
//             self.f_n_x0[1]

void pCoeff(std::vector<RCP<const Basic>> &arr, std::size_t order) {
  std::vector<std::vector<RCP<const Basic>>> dict(
      order + 1, std::vector<RCP<const Basic>>(order + 1));
  for (std::size_t j{1}; j <= order; ++j) {
    dict[j][j] = pow(arr[1], SymEngine::integer(j));
  }
  printDict(dict);
  // std::cout << dict.size();
}

int main() {
  Timer t;
  SymEngine::print_stack_on_segfault();

  RCP<const Symbol> x = symbol("x");
  RCP<const Symbol> rs = symbol("ra");
  auto ex = add(add(x, rs), log(sub(add(x, rs), integer(1))));

  std::vector<RCP<const Basic>> ser{fourierSeries(ex, x, 3)};
  // printArray(ser);
  pCoeff(ser, 3);
  std::cout << "Time elapsed: " << t.elapsed() << " seconds\n";
  return 0;
}
