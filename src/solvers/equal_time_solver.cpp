#include <future>
#include <boost/math/tools/roots.hpp>
#include "equal_time_solver.hpp"

using boost::math::tools::toms748_solve;
using boost::math::tools::eps_tolerance;

EqualTimeSolver& EqualTimeSolver::withRestriction(EnergyDensity& rho3)
{
    restriction = rho3;
    hasRestriction = true;
    return *this; 
}

HighPrecision EqualTimeSolver::getEqualTime()
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
    
        const int digits = std::numeric_limits<HighPrecision>::digits;
        HighPrecision fa = h(lowerLimit);
        HighPrecision fb = h(upperLimit);
    
        if (fa * fb >= 0) {
            throw std::runtime_error("No root found.");
        }
    
        std::uintmax_t max_iter = 10;
        auto result = toms748_solve(h, lowerLimit, upperLimit, fa, fb, eps_tolerance<HighPrecision>(digits), max_iter);
        HighPrecision timeEquality = (result.first + result.second) / HighPrecision(2.0);
    
        //HighPrecision rho1Equal = rho1(p, timeEquality);
        HighPrecision rho2Equal = rho2(timeEquality);
        // Check for constraint
        if (hasRestriction == true && !(rho2Equal > restriction(timeEquality)))
        {
            throw std::runtime_error("No time of equality with the given constraint.");
        }
        return timeEquality;
}