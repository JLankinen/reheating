#include <future>
#include <boost/math/tools/roots.hpp>
#include "solvers/equal_time_solver.hpp"
#include <iostream>
#include <iomanip>
#include <algorithm>

using boost::math::tools::toms748_solve;
using boost::math::tools::eps_tolerance;

struct Bracket
{
    HighPrecision low;
    HighPrecision high;
    HighPrecision fa;
    HighPrecision fb;
};

/**
 * Custom bracketing function. Increases upper limit by times 10^1 until sign change is found (this will
 * happen at some point) and returns the found bracket to be used in Toms method.
 */
Bracket findBracket(EnergyDensity h, HighPrecision low)
{
    HighPrecision fa = h(low);
    HighPrecision high = low * HighPrecision("10");
    HighPrecision fb = h(high);

    const int maxAttempts = 150;
    int attempts = 0;
    // Increase upperlimit until bracket is found. Equations show this will happen at some point.
    
    while (fa * fb > 0 && attempts < maxAttempts)
    {
        high *= HighPrecision("10");
        fb = h(high);
        attempts++;
    }

    if (fa * fb > 0)
    {
        throw std::runtime_error("Failed to bracket root.");
    }

    
    return {low, high, fa, fb};
}

std::tuple<HighPrecision, HighPrecision, HighPrecision> EqualTimeSolver::getEqualTime()
{
        // Function difference

        auto h = [&](HighPrecision t) -> HighPrecision {
            HighPrecision val1 = rho1(t);
            HighPrecision val2 = rho2(t);

            if (val1 <= 0 || val2 <= 0)
            {
                return val1 - val2;
            }

            HighPrecision log1 = log(val1);
            HighPrecision log2 = log(val2);
            HighPrecision result = log1 - log2;
            //std::cout << "  [TOMS748] Evaluating h(t): t = " << t << ", log_rho1(t) = " << log1 << ", log_rho2(t)= " << log2 <<", h(t)= " << result << "\n";
            return result;
        };


        std::pair<HighPrecision, HighPrecision> result;
        const int digits = std::numeric_limits<HighPrecision>::digits;
        Bracket bracket = findBracket(h, lowerLimit);
        result = toms748_solve(h, bracket.low, bracket.high, bracket.fa, bracket.fb, eps_tolerance<HighPrecision>(digits), maxIter);

        HighPrecision timeEquality = (result.first + result.second) / HighPrecision("2.0");
        HighPrecision rho1Equal = rho1(timeEquality);
        HighPrecision rho2Equal = rho2(timeEquality);

        return std::make_tuple(timeEquality, rho1Equal, rho2Equal);
}


// New logic

std::optional<std::tuple<HighPrecision, HighPrecision, HighPrecision>> EqualTimeSolver::findEqualTime()
{
    const int digits = std::numeric_limits<HighPrecision>::digits;
    auto h = [&](HighPrecision t) -> HighPrecision
    {
        HighPrecision val1 = rho1(t);
        HighPrecision val2 = rho2(t);

        if (val1 <= 0 | val2 <= 0)
        {
            return val1 - val2;
        }

        HighPrecision log1 = log(val1);
        HighPrecision log2 = log(val2);
        HighPrecision result = log1 - log2;
        //std::cout << "  [TOMS748] Evaluating h(t): t = " << t << ", log_rho1(t) = " << log1 << ", log_rho2(t)= " << log2 <<", h(t)= " << result << "\n";
        return result;
    };

    Bracket bracket = findBracket(h, lowerLimit);

    if (bracket.fa * bracket.fb > 0)
    {
        return std::nullopt; // No root exists
    }

    try
    {
        auto result = toms748_solve(h, bracket.low, bracket.high, bracket.fa, bracket.fb, eps_tolerance<HighPrecision>(digits), maxIter);
        HighPrecision timeEquality = (result.first + result.second) / HighPrecision("2.0");
        HighPrecision rho1Equal = rho1(timeEquality);
        HighPrecision rho2Equal = rho2(timeEquality);

        return std::make_tuple(timeEquality, rho1Equal, rho2Equal);
    }
    catch (...)
    {
        return std::nullopt;
    }

}