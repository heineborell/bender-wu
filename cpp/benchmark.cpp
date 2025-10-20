#include <chrono>
#include <iostream>
#include <symengine/derivative.h>
#include <symengine/functions.h>
#include <symengine/integer.h>
#include <symengine/pow.h>
#include <symengine/symbol.h>

using namespace SymEngine;

int main() {
  auto r = symbol("r");

  // =============================
  // Case 1: Abstract symbolic V(r)
  // =============================
  vec_basic args{r}; // FIX: explicitly use vec_basic
  auto V = function_symbol("V", args);

  auto t1 = std::chrono::high_resolution_clock::now();
  RCP<const Basic> dV = V;
  for (int i = 0; i < 50; i++) {
    dV = diff(dV, r); // Derivative(V(r), r, r, ...)
  }
  auto t2 = std::chrono::high_resolution_clock::now();

  std::cout << "Case 1 result: " << *dV << std::endl;
  std::cout << "Case 1 time: "
            << std::chrono::duration<double, std::milli>(t2 - t1).count()
            << " ms\n";

  // =============================
  // Case 2: Explicit function r^100
  // =============================
  auto f = pow(r, integer(100));

  auto t3 = std::chrono::high_resolution_clock::now();
  RCP<const Basic> df = f;
  for (int i = 0; i < 50; i++) {
    df = diff(df, r); // compute actual derivative each time
  }
  auto t4 = std::chrono::high_resolution_clock::now();

  std::cout << "Case 2 result: " << *df << std::endl;
  std::cout << "Case 2 time: "
            << std::chrono::duration<double, std::milli>(t4 - t3).count()
            << " ms\n";

  return 0;
}
