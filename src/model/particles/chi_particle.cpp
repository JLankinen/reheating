#include <cmath>
#include <boost/math/quadrature/gauss_kronrod.hpp>

#include "chi_particle.hpp"
#include "model/energy/creation_decay.hpp"

using boost::math::quadrature::gauss_kronrod;

EnergyDensity ChiParticle::energyDensityStiff()
{
    return [this](HighPrecision t) -> HighPrecision
    {
        constexpr double n = 1.0;  // For stiff matter universe, n = 1.
        constexpr double tol = 1e-15;
        
        HighPrecision prefactor = 1 / pow(t, HighPrecision(4.0 / 3.0));
        
        auto integrand = [&](HighPrecision tprime)
        {
            HighPrecision val = ChiDecayRate(p, n, p.t0, tprime) * rhoPhi(tprime)
                                * pow(t, HighPrecision(4.0) / HighPrecision(3.0));
            return val;
        };
    
        static gauss_kronrod<HighPrecision, 10> integrator;
        HighPrecision integralResult = integrator.integrate(integrand, p.t0, t, tol);
        return prefactor * integralResult;
    };
}

