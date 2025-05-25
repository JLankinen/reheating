#include <map>
#include <mutex>
#include <boost/math/quadrature/gauss_kronrod.hpp>

#include "phi_particle.hpp"
#include "utils/types.hpp"

#include "model/energy/creation_decay.hpp"

using boost::math::quadrature::gauss_kronrod;

std::map<HighPrecision, HighPrecision> rhoPhiCache;
std::mutex cacheMutex;

HighPrecision quantize(const HighPrecision& value, int digits = 30) {
    std::ostringstream oss;
    oss << std::setprecision(digits) << std::fixed << value;
    return HighPrecision(oss.str());
}


EnergyDensity PhiParticle::energyDensityStiff()
{
    return [this](HighPrecision t) -> HighPrecision
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
    };
}