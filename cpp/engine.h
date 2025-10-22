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

std::vector<RCP<const Basic>>
fourierSeries(RCP<const Basic> func, RCP<const Symbol> var, std::size_t order);

std::vector<std::vector<RCP<const Basic>>>
inverseCoeff(std::vector<RCP<const Basic>> &arr, std::size_t order);

RCP<const Basic> getPotential(RCP<const Symbol> r, const int L, const int s);
