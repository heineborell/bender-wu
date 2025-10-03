#include <chrono>

#include <cstddef>
#include <symengine/add.h>
#include <symengine/basic.h>
#include <symengine/functions.h>
#include <symengine/pow.h>
#include <symengine/rational.h>
#include <symengine/series.h>
#include <symengine/symengine_casts.h>
#include <symengine/symengine_rcp.h>

using SymEngine::add;
using SymEngine::Basic;
using SymEngine::integer;
using SymEngine::log;
using SymEngine::RCP;
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

void fourierSeries(std::vector<RCP<const Basic>> &arr, RCP<const Basic> func,
                   RCP<const Symbol> var, std::size_t order) {
  arr.push_back(func);            // 0.th order
  arr.push_back(func->diff(var)); // 1. order
  for (std::size_t i{2}; i <= order; ++i) {
    arr.push_back(arr[i - 1]->diff(var));
  }
}

void printArray(std::vector<RCP<const Basic>> &arr) {
  for (RCP<const Basic> item : arr) {
    std::cout << *item << '\n';
  }
}

int main() {
  Timer t;
  SymEngine::print_stack_on_segfault();

  RCP<const Symbol> x = symbol("x");
  RCP<const Symbol> a = symbol("a");
  std::vector<RCP<const Basic>> arr{};
  auto ex = add(add(x, a), log(sub(add(x, a), integer(1))));

  fourierSeries(arr, ex, x, 30);
  printArray(arr);
  std::cout << "Time elapsed: " << t.elapsed() << " seconds\n";
  return 0;
}
