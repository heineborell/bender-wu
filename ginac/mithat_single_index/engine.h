#pragma once

// Disable *all* warnings for GiNaC on Clang (macOS) since macOS is clang is
// more strict than the Linux g++
#if defined(__clang__)
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Weverything"
#endif

#ifdef IN_GINAC
#include "ginac.h"
#else
#include <ginac/ginac.h>
#endif

#if defined(__clang__)
#pragma clang diagnostic pop
#endif

#include <vector>
using namespace GiNaC;

inline constexpr int nu{3};
inline constexpr int expansionOrder{200};
inline constexpr int lmax{2 * expansionOrder};
inline int even{0};
inline int lstep{2};

std::vector<ex> fourierSeries(const ex &func, const symbol &var, int order);
std::vector<ex> vSeries(const ex &func, const symbol &var, int order);
std::vector<std::vector<ex>> pCoeff(const std::vector<ex> &arr,
                                    std::size_t order);
std::vector<ex> cCoeff(const std::vector<ex> &f_n,
                       const std::vector<std::vector<ex>> &dict, const ex &func,
                       const symbol &var, std::size_t order);
void evenCheck(const std::vector<ex> &f_n_x0);
std::vector<ex> energy(const std::vector<ex> &f_n_x0, ex &omega);
