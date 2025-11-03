#ifdef IN_GINAC
#include "ginac.h"
#else
#include <ginac/ginac.h>
#endif

using namespace GiNaC;

std::vector<ex> fourierSeries(const ex &func, symbol &var, std::size_t order);
std::vector<std::vector<ex>> pCoeff(std::vector<ex> &arr, std::size_t order);
