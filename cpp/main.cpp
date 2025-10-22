#include "engine.h"
#include "helper.h"
#include <cstddef>

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

  // std::cout << "solution is " << *solve(potential->diff(r), r) << '\n';
  // auto solutions = solve(potential->diff(r), r);

  std::vector<RCP<const Basic>> ser{fourierSeries(ex, x, order + 1)};
  printArray(ser);

  std::vector<std::vector<RCP<const Basic>>> inverse_coeff{
      inverseCoeff(ser, order)};
  // printDict(inverse_coeff);
  // auto contain = solutions.get_container();
  // auto finite_set =
  //     SymEngine::rcp_static_cast<const SymEngine::FiniteSet>(solutions);
  // std::cout << "the set is " << *finite_set;
  // auto contain = finite_set->get_container();
  // for (auto &set : finite_set) {
  //   std::cout << *set << '\n';
  // }
  std::cout << "Time elapsed: " << t.elapsed() << " seconds\n";
  return 0;
}
