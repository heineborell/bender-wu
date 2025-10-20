#include <chrono>

#include <cstddef>
#include <symengine/add.h>
#include <symengine/basic.h>
#include <symengine/functions.h>
#include <symengine/integer.h>
#include <symengine/ntheory.h>
#include <symengine/pow.h>
#include <symengine/rational.h>
#include <symengine/series.h>
#include <symengine/sets.h>
#include <symengine/simplify.h>
#include <symengine/solve.h>
#include <symengine/symengine_casts.h>
#include <symengine/symengine_rcp.h>

using SymEngine::add;
using SymEngine::Basic;
using SymEngine::FiniteSet;
using SymEngine::finiteset;
using SymEngine::integer;
using SymEngine::log;
using SymEngine::max;
using SymEngine::Max;
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
      if (dict[i][j].get()) // check if its nullptr (since these are smart gotta
                            // use .get)
        std::cout << '(' << i << ',' << j << "): " << *dict[i][j] << '\n';
      else
        std::cout << '(' << i << ',' << j << "): " << "nullptr" << '\n';
    }
  }
}
void printArray(std::vector<RCP<const Basic>> &arr) {
  for (const RCP<const Basic> &item : arr) {
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

int main() {
  Timer t;
  SymEngine::print_stack_on_segfault();

  RCP<const Symbol> x = symbol("x");
  RCP<const Symbol> r = symbol("r");
  RCP<const Symbol> rs = symbol("ra");

  const int order{3};
  const int N{3};
  const int s{2};
  const int L{2};

  auto ex = add(add(x, rs), log(sub(add(x, rs), integer(1))));
  auto potential = getPotential(r, L, s);

  std::cout << "solution is " << *solve(potential->diff(r), r) << '\n';
  auto solutions = solve(potential->diff(r), r);

  std::vector<RCP<const Basic>> ser{fourierSeries(ex, x, order + 1)};
  std::vector<std::vector<RCP<const Basic>>> inverse_coeff{
      inverseCoeff(ser, order)};
  printDict(inverse_coeff);
  // auto contain = solutions.get_container();
  auto finite_set =
      SymEngine::rcp_static_cast<const SymEngine::FiniteSet>(solutions);
  std::cout << "the set is " << *finite_set;
  // auto contain = finite_set->get_container();
  // for (auto &set : finite_set) {
  //   std::cout << *set << '\n';
  // }
  std::cout << "Time elapsed: " << t.elapsed() << " seconds\n";
  return 0;
}
