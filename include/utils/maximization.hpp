#ifndef MAXIMIZATION_H_
#define MAXIMIZATION_H_

#include <boost/math/tools/minima.hpp>
#include "utils/types.hpp"

using boost::math::tools::brent_find_minima;

template<typename F>
double maximize(F&& f, double min, double max)
{
    const int digits = std::numeric_limits<double>::digits;
    auto result = brent_find_minima(
        [&](double t) {return -(std::forward<F>(f))(t);},
        min,
        max,
        digits);
    double location = result.first;
    return location;
}

#endif

