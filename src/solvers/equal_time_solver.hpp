#ifndef EQUAL_TIME_SOLVER_H_
#define EQUAL_TIME_SOLVER_H_

#include <functional>
#include "parameters/parameters.hpp"

using EnergyDensity = std::function<HighPrecision(ModelParameters, HighPrecision t)>;

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
    private:
        EnergyDensity rho1;
        EnergyDensity rho2;
        EnergyDensity restriction;
        ModelParameters p;
        HighPrecision lowerLimit;
        HighPrecision upperLimit;
        bool hasRestriction;
    public:
        EqualTimeSolver(EnergyDensity _rho1, EnergyDensity _rho2, ModelParameters _p, HighPrecision _lowerLimit, HighPrecision _upperLimit):
        rho1{_rho1}, rho2{_rho2}, p{_p}, lowerLimit{_lowerLimit}, upperLimit{_upperLimit}, hasRestriction{false} {};

        EqualTimeSolver& withRestriction(EnergyDensity rho3);
        HighPrecision getEqualTime();

};

#endif