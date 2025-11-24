#pragma once
#ifdef IN_GINAC
#include "ginac.h"
#else
#include <ginac/ginac.h>
#endif

using namespace GiNaC;

inline constexpr int nu{2};
inline constexpr int eLL{10};

std::vector<ex> fourierSeries(const ex &func, const symbol &var, int order);
std::vector<std::vector<ex>> pCoeff(const std::vector<ex> &arr,
                                    std::size_t order);
std::vector<ex> cCoeff(const std::vector<ex> &f_n,
                       const std::vector<std::vector<ex>> &dict, const ex &func,
                       const symbol &var, std::size_t order);
std::vector<ex> aCoeff(const std::vector<ex> &f_n_x0, ex &omega, int order);
std::vector<ex> Energy(std::vector<std::vector<ex>> &A, std::vector<ex> &f_n_x0,
                       ex &omega);
