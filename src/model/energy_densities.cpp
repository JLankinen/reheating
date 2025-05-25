#include <boost/multiprecision/cpp_dec_float.hpp>
#include <cmath>
#include <numbers>
#include <functional>
#include <map>
#include <mutex>
#include <future>
#include <iostream>
#include <boost/math/quadrature/gauss_kronrod.hpp>
#include <boost/math/tools/roots.hpp>
#include "energy_densities.hpp"
#include "parameters/parameters.hpp"
#include "model/creation_decay.hpp"

using boost::math::quadrature::gauss_kronrod;
using boost::math::tools::toms748_solve;
using boost::math::tools::eps_tolerance;

using HighPrecision = boost::multiprecision::cpp_dec_float_100;

std::map<HighPrecision, HighPrecision> rhoPhiCache;
std::mutex cacheMutex;

HighPrecision quantize(const HighPrecision& value, int digits = 30) {
    std::ostringstream oss;
    oss << std::setprecision(digits) << std::fixed << value;
    return HighPrecision(oss.str());
}


HighPrecision RhoStiff(ModelParameters p, HighPrecision t)
{
    return HighPrecision(1) / (HighPrecision(24.0) * std::numbers::pi * p.G_N * pow(t, HighPrecision(2.0)));
}


HighPrecision RhoPhiStiff(ModelParameters p, HighPrecision t)
{
    constexpr double n = 1.0;
    constexpr double tol = 1e-15;

    HighPrecision key = quantize(t, 30);  // Quantize to 30 digits for cache key

    {
        std::lock_guard<std::mutex> lock(cacheMutex);
        auto it = rhoPhiCache.find(key);
        if (it != rhoPhiCache.end()) {
            return it->second;
        }
    }
    
    HighPrecision prefactor = (HighPrecision(1.0) / t) * exp(-ChiDecayRate(p, n, p.t0, t));

    auto integrand = [&](HighPrecision tprime)
    {
        HighPrecision creationRate = PhiCreationRate(p, tprime);
        HighPrecision decayRate = ChiDecayRate(p, n, p.t0, tprime);
        HighPrecision logVal = log(tprime) + log(creationRate) + decayRate;
        HighPrecision val = exp(logVal);
        return val;
    };

    static gauss_kronrod<HighPrecision, 10> integrator;
    HighPrecision integralResult = integrator.integrate(integrand, p.t0, t, tol);
    HighPrecision result = prefactor * integralResult;

    {
        std::lock_guard<std::mutex> lock(cacheMutex);
        rhoPhiCache[key] = result;
    }

    return result;
}


HighPrecision RhoChiStiff(ModelParameters p, HighPrecision t)
{
    constexpr double n = 1.0;  // For stiff matter universe, n = 1.
    constexpr double tol = 1e-15;
    
    HighPrecision prefactor = 1 / pow(t, HighPrecision(4.0 / 3.0));
    
    auto integrand = [&](HighPrecision tprime)
    {
        HighPrecision val = ChiDecayRate(p, n, p.t0, tprime) * RhoPhiStiff(p, tprime) * pow(t, HighPrecision(4.0) / HighPrecision(3.0));
        return val;
    };

    static gauss_kronrod<HighPrecision, 10> integrator;
    HighPrecision integralResult = integrator.integrate(integrand, p.t0, t, tol);
    return prefactor * integralResult;
}


HighPrecision EqualTime(std::function<HighPrecision(ModelParameters, HighPrecision t)> f1,
                 std::function<HighPrecision(ModelParameters, HighPrecision t)> f2,
                 std::function<HighPrecision(ModelParameters, HighPrecision t)> restrictionFunc,
                 ModelParameters p,
                 HighPrecision lowerLimit,
                 HighPrecision upperLimit)
{

    // Function difference
    auto h = [&](HighPrecision t) -> HighPrecision {
        // Run f1 and f2 asynchronously
        auto val1_fut = std::async(std::launch::async, f1, p, t);
        auto val2_fut = std::async(std::launch::async, f2, p, t);
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

    // Check for constraint
    if(f2(p, timeEquality) > restrictionFunc(p, timeEquality))
    {
        return timeEquality;
    }
    else
    {
        throw std::runtime_error("No time of equality with the given constraint.");
    }
}
