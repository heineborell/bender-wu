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
inline constexpr int expansionOrder{50};
inline constexpr int lmax{2 * expansionOrder};
inline int even{0};
inline int lstep{2};

std::array<ex, 2 * expansionOrder + 1> fourierSeries(const ex &func,
                                                     const symbol &var);
std::array<ex, 2 * expansionOrder + 1> vSeries(const ex &func,
                                               const symbol &var);
std::vector<std::vector<ex>> pCoeff(const std::vector<ex> &arr,
                                    std::size_t order);
std::vector<ex> cCoeff(const std::vector<ex> &f_n,
                       const std::vector<std::vector<ex>> &dict, const ex &func,
                       const symbol &var, std::size_t order);
void evenCheck(const std::vector<ex> &f_n_x0);
std::array<ex, lmax / 2 + 1>
energy(const std::array<ex, 2 * expansionOrder + 1> &f_n_x0, ex &omega);
