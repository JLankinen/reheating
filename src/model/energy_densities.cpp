#include <boost/multiprecision/cpp_dec_float.hpp>
#include <cmath>
#include <numbers>
#include <functional>
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

ModelParameters p;

template <typename F, typename ...Args>
auto TimeCounterDecorator(F&& f, Args&&... args)
{
    auto start = std::chrono::high_resolution_clock::now();
    auto result = std::invoke(std::forward<F>(f), std::forward<Args>(args)...);
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    std::cout << "Elapsed time: " << duration.count()  << " µs" << std::endl;
    return result;
}

HighPrecision RhoStiff(HighPrecision t)
{
    return HighPrecision(1) / (HighPrecision(24.0) * std::numbers::pi * p.G_N * pow(t, HighPrecision(2.0)));
}


HighPrecision RhoPhiStiff(ModelParameters p, HighPrecision t)
{
    constexpr double n = 1.0;  // For stiff matter universe, n = 1.
    constexpr double tol = 1e-8;
    
    HighPrecision prefactor = (HighPrecision(1.0) / t) * exp(-ChiDecayRate(p, n, p.t0, t));

    auto integrand = [&](HighPrecision tprime)
    {
        HighPrecision creationRate = PhiCreationRate(p, tprime);
        HighPrecision decayRate = ChiDecayRate(p, n, p.t0, tprime);
        HighPrecision val = tprime * creationRate * exp(decayRate);
        std::cout << "RhoPhiStiff:\t" << "tprime: " << tprime << ", integrand: " << val << std::endl;
        return val;
    };

    static gauss_kronrod<HighPrecision, 31> integrator;
    HighPrecision integralResult = integrator.integrate(integrand, p.t0, t, tol);
    std::cout << "RhoPhiStiff Integral: " << prefactor * integralResult << std::endl;
    return prefactor * integralResult;

}


HighPrecision RhoChiStiff(ModelParameters p, HighPrecision t)
{
    constexpr double n = 1.0;  // For stiff matter universe, n = 1.
    constexpr double tol = 1e-8;
    
    HighPrecision prefactor = 1 / pow(t, HighPrecision(4.0 / 3.0));
    
    auto integrand = [&](HighPrecision tprime)
    {
        return ChiDecayRate(p, n, p.t0, tprime) * RhoPhiStiff(p, tprime) * pow(t, HighPrecision(4.0) / HighPrecision(3.0));
    };

    static gauss_kronrod<HighPrecision, 31> integrator;
    HighPrecision integralResult = integrator.integrate(integrand, p.t0, t, tol);
    return prefactor * integralResult;
}


HighPrecision EqualTime(std::function<HighPrecision(HighPrecision t)> f1,
                 std::function<HighPrecision(ModelParameters, HighPrecision t)> f2,
                 std::function<HighPrecision(ModelParameters, HighPrecision t)> restrictionFunc,
                 ModelParameters p,
                 HighPrecision lowerLimit,
                 HighPrecision upperLimit)
{

    // Function difference
    auto h = [&](HighPrecision t) -> HighPrecision {
        return f1(t) - f2(p, t);
    };

    const int digits = std::numeric_limits<HighPrecision>::digits;
    HighPrecision fa = h(lowerLimit);
    HighPrecision fb = h(upperLimit);
    std::uintmax_t max_iter = 10;
    auto result = toms748_solve(h, lowerLimit, upperLimit, fa, fb, eps_tolerance<HighPrecision>(digits), max_iter);
    auto timeEquality = (result.first + result.second) / 2.0;

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
