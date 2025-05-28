#include <future>
#include <boost/math/tools/roots.hpp>
#include "solvers/equal_time_solver.hpp"
#include <iostream>
#include <iomanip>
#include <algorithm>

using boost::math::tools::toms748_solve;
using boost::math::tools::eps_tolerance;

EqualTimeSolver& EqualTimeSolver::maxIterations(std::uintmax_t max)
{
    maxIter = max;
    return *this;
}


EqualTimeSolver& EqualTimeSolver::withRestriction(EnergyDensity rho3)
{
    restriction = rho3;
    hasRestriction = true;
    return *this; 
}


std::tuple<HighPrecision, HighPrecision, HighPrecision> EqualTimeSolver::getEqualTime()
{
        // Function difference
        auto h = [&](HighPrecision t) -> HighPrecision {
            // Run rho1 and rho2 asynchronously
            auto val1_fut = std::async(std::launch::async, rho1, t);
            auto val2_fut = std::async(std::launch::async, rho2, t);
            HighPrecision val1 = val1_fut.get();
            HighPrecision val2 = val2_fut.get();
            return val1 - val2;
        };
        auto [low, high] = findBracket(h, lowerLimit);
        const int digits = std::numeric_limits<HighPrecision>::digits;
        HighPrecision fa = h(low);
        HighPrecision fb = h(high);

        auto result = toms748_solve(h, low, high, fa, fb, eps_tolerance<HighPrecision>(digits), maxIter);
        HighPrecision timeEquality = (result.first + result.second) / HighPrecision(2.0);
    
        HighPrecision rho1Equal = rho1(timeEquality);
        HighPrecision rho2Equal = rho2(timeEquality);
        // Check for constraint
        if (hasRestriction == true && !(rho2Equal > restriction(timeEquality)))
        {
            throw std::runtime_error("No time of equality with the given constraint.");
        }
        return std::make_tuple(timeEquality, rho1Equal, rho2Equal);
}


std::pair<HighPrecision, HighPrecision> EqualTimeSolver::findBracket(const EnergyDensity& rho,
                                                                     HighPrecision t0,
                                                                     HighPrecision initialStep,
                                                                     HighPrecision maxStep,
                                                                     HighPrecision growthFactor,
                                                                     int maxAttempts)
{
    HighPrecision a = t0;
    HighPrecision fa = rho(a);
    HighPrecision step = initialStep;
    HighPrecision b = a + step;
    HighPrecision fb = rho(b);
    int attempts = 0;

    // Find the bracket by increasing the interval in steps.
    while(fa * fb > 0 && attempts <= maxAttempts && step <= maxStep)
    {
        // Change bracket location
        a = b;
        fa = fb;
        // Increase the right side.
        step *= growthFactor;
        b = a + step;
        fb = rho(b);
        ++attempts;
    }

    if (fa * fb > 0)
    {
        throw std::runtime_error("Failed to bracket root.");
    }

    return {a, b};
}

