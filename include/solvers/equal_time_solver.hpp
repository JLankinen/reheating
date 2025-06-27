#ifndef EQUAL_TIME_SOLVER_H_
#define EQUAL_TIME_SOLVER_H_

#include <tuple>
#include <utility>
#include "parameters/parameters.hpp"
#include "utils/types.hpp"
/**
 * @class EqualTimeSolver
 * @brief Finds the time when two energy density functions are equal within a specified interval.
 *
 * This class encapsulates the root-finding logic to determine the time `t` in a given interval
 * [lowerLimit, upperLimit] where two provided energy density functions (rho1 and rho2) are equal,
 * i.e., rho1(t) == rho2(t). The root-finding uses the Toms748 algorithm from Boost.Math.
 *
 * Optionally, a restriction function can be supplied to impose a constraint on the solution,
 * such as requiring that rho2(t) exceeds a threshold function at the root.
 *
 * The solver performs asynchronous evaluation of rho1 and rho2 for improved performance.
 */

class EqualTimeSolver
{
    public:
        enum class Method
        {
            Toms748,
            Bisect
        };
    private:
        EnergyDensity rho1;
        EnergyDensity rho2;
        HighPrecision lowerLimit;
        std::uintmax_t maxIter = 20;
    public:
        EqualTimeSolver(EnergyDensity _rho1, EnergyDensity _rho2, HighPrecision _lowerLimit):
        rho1{_rho1}, rho2{_rho2}, lowerLimit{_lowerLimit} {};
       
        /**
         * @brief Get the time t_eq when rho1 and rho2 are equal i.e., rho1(t)=rho2(t).
         * 
         * @return A tuple consisting of (t_eq, rho1(t_eq), rho2(t_eq))
         */
        std::tuple<HighPrecision, HighPrecision, HighPrecision> getEqualTime();

}; 

#endif