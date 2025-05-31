#ifndef MAXIMIZATION_H_
#define MAXIMIZATION_H_

#include <boost/math/tools/minima.hpp>
#include "utils/types.hpp"

using boost::math::tools::brent_find_minima;

template<typename F>
HighPrecision maximize(F&& f, HighPrecision min, HighPrecision max)
{
    const int digits = std::numeric_limits<HighPrecision>::digits;
    auto result = brent_find_minima(
        [&](HighPrecision t) {return -(std::forward<F>(f))(t);},
        min,
        max,
        digits);
    HighPrecision location = result.first;
    return location;
}

#endif

