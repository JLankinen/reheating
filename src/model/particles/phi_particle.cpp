#include <cmath>

#include "model/particles/phi_particle.hpp"
#include "model/energy/creation_decay.hpp"
#include "utils/types.hpp"
#include "utils/integration.hpp"


std::map<HighPrecision, HighPrecision> PhiParticle::rhoPhiCache;
std::mutex PhiParticle::cacheMutex;

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
        
        HighPrecision prefactor = (HighPrecision(1.0) / t) * exp(-ChiDecayRate(this->p, n, this->p.t0, t));

        auto integrand = [&](HighPrecision tprime)
        {
            HighPrecision creationRate = PhiCreationRate(this->p, tprime);
            HighPrecision decayRate = ChiDecayRate(this->p, n, this->p.t0, tprime);
            HighPrecision logVal = log(tprime) + log(creationRate) + decayRate;
            HighPrecision val = exp(logVal);
            return val;
        };

        HighPrecision integralResult = IntegrationUtils::integrate(integrand, this->p.t0, t, tol);
        HighPrecision result = prefactor * integralResult;

        {
            std::lock_guard<std::mutex> lock(cacheMutex);
            rhoPhiCache[key] = result;
        }

        return result;
    };
}


EnergyDensity PhiParticle::energyDensityMatter(HighPrecision t0, HighPrecision initialRho)
{
    return [this, t0, initialRho](HighPrecision t)->HighPrecision{
        // n = 4 in matter dominated Universe
        constexpr double n = 4.0;
        HighPrecision time = pow(t0, HighPrecision(2.0 / 3.0)) / pow(t, HighPrecision(2.0 / 3.0));
        return pow(time, HighPrecision(3)) * exp(-ChiDecayRate(this->p, n, t0, t)) * initialRho;
    };
}


EnergyDensity PhiParticle::energyDensityRadiation(HighPrecision t0, HighPrecision initialRho)
{
    return [this, t0, initialRho](HighPrecision t)->HighPrecision{
        // n = 2 in matter dominated Universe
        constexpr double n = 2.0;
        HighPrecision time = pow(t0, HighPrecision(1.0 / 2.0)) / pow(t, HighPrecision(1.0 / 2.0));
        return pow(time, HighPrecision(3)) * exp(-ChiDecayRate(this->p, n, t0, t)) * initialRho;
    };    
}