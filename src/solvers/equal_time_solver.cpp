#include <future>
#include <boost/math/tools/roots.hpp>
#include "solvers/equal_time_solver.hpp"
#include <iostream>
#include <iomanip>
#include <algorithm>

using boost::math::tools::toms748_solve;
using boost::math::tools::bisect;
using boost::math::tools::eps_tolerance;


EqualTimeSolver& EqualTimeSolver::usingMethod(Method meth)
{
    method = meth;
    return *this;
}

void EqualTimeSolver::setUpperLimit(HighPrecision ul)
{
    this->upperLimit = ul;
}


std::tuple<HighPrecision, HighPrecision, HighPrecision> EqualTimeSolver::getEqualTime()
{
    switch (method)
    {
    case Method::Toms748:
        return toms748Method();
    case Method::Bisect:
        return bisectMethod();
    default:
        throw std::runtime_error("No root finding method selected.");
    }
}

/**
 * Custom bracketing function. Increases upper limit by times 10^1 until sign change is found (this will
 * happen at some point) and returns the found bracket to be used in Toms method.
 */
std::tuple<HighPrecision, HighPrecision> findBracket(EnergyDensity h, HighPrecision low)
{
    HighPrecision fa = h(low);
    HighPrecision high = low * HighPrecision("1e3");
    HighPrecision fb = h(high);
    const int maxAttempts = 100;
    int attempts = 0;
    // Increase upperlimit until bracket is found. Equations show this will happen at some point.
    if (fa > 0)
    {
        while (fb > 0 && attempts < maxAttempts)
        {
            high = high * HighPrecision("1e2");
            fb = h(high);
            attempts++;
        }
    }
    else if (fa < 0)
    {
        while (fb < 0 && attempts < maxAttempts)
        {
            high = high * HighPrecision("1e2");
            fb = h(high);
            attempts++;
        }
    }
    return std::make_tuple(low, high);
}

std::tuple<HighPrecision, HighPrecision, HighPrecision> EqualTimeSolver::toms748Method()
{
        // Function difference
        auto h = [&](HighPrecision t) -> HighPrecision {
            HighPrecision val1 = rho1(t);
            HighPrecision val2 = rho2(t);
            HighPrecision result = val1 - val2;
            std::cout << "  [TOMS748] Evaluating h(t): t = " << t << ", rho1(t) = " << val1 << ", rho2(t)= " << val2 <<", h(t)= " << result << "\n";
            return result;
        };


        std::pair<HighPrecision, HighPrecision> result;
        //auto low = lowerLimit;
        //auto high = upperLimit;
        const int digits = std::numeric_limits<HighPrecision>::digits;
        auto [low, high] = findBracket(h, lowerLimit);
        HighPrecision fa = h(low);
        HighPrecision fb = h(high);
        result = toms748_solve(h, low, high, fa, fb, eps_tolerance<HighPrecision>(digits), maxIter);

        HighPrecision timeEquality = (result.first + result.second) / HighPrecision(2.0);
        HighPrecision rho1Equal = rho1(timeEquality);
        HighPrecision rho2Equal = rho2(timeEquality);

        return std::make_tuple(timeEquality, rho1Equal, rho2Equal);
}


std::tuple<HighPrecision, HighPrecision, HighPrecision> EqualTimeSolver::bisectMethod()
{
    // Function difference
    auto h = [&](HighPrecision t) -> HighPrecision {
        HighPrecision val1 = rho1(t);
        HighPrecision val2 = rho2(t);
        HighPrecision result = val1 - val2;
        std::cout << "  [Bisect] Evaluating h(t): t = " << t << ", rho1(t) = " << val1 << ", rho2(t)= " << val2 <<", h(t)= " << result << "\n";
        return result;
    };

    const int digits = std::numeric_limits<HighPrecision>::digits;
    std::uintmax_t maxIter = 100;
    auto result = bisect(h, lowerLimit, upperLimit, eps_tolerance<HighPrecision>(digits), maxIter);
    HighPrecision eqTime = (result.first + result.second) / 2;

    HighPrecision rho1Equal = rho1(eqTime);
    HighPrecision rho2Equal = rho2(eqTime);

    return std::make_tuple(eqTime, rho1Equal, rho2Equal);
}