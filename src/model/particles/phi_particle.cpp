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

void PhiParticle::setInitialRhoMatter(const HighPrecision& rhoInit)
{
    this->initialRhoMatter = rhoInit;
}

HighPrecision PhiParticle::getInitialRhoMatter() const
{
    return initialRhoMatter;
}

void PhiParticle::setInitialRhoRadiation(const HighPrecision& rhoInit)
{
    this->initialRhoRadiation = rhoInit;
}

HighPrecision PhiParticle::getInitialRhoRadiation() const
{
    return initialRhoRadiation;
}



EnergyDensity PhiParticle::energyDensityStiff()
{
    constexpr double n = 1.0;
    ChiDecayRate chiDecay(this->p, n, this->p.t0);

    return [this, chiDecay](HighPrecision t) -> HighPrecision
    {

        HighPrecision key = quantize(t, 30);  // Quantize to 30 digits for cache key

        {
            std::lock_guard<std::mutex> lock(cacheMutex);
            auto it = rhoPhiCache.find(key);
            if (it != rhoPhiCache.end()) {
                return it->second;
            }
        }
        
        HighPrecision prefactor = (HighPrecision(1.0) / t) * exp(-chiDecay(t));

        auto integrand = [&](HighPrecision tprime)
        {
            HighPrecision creationRate = PhiCreationRate(this->p, tprime);
            HighPrecision decayRate = chiDecay(tprime);
            HighPrecision logVal = log(tprime) + log(creationRate) + decayRate;
            HighPrecision val = exp(logVal);
            return val;
        };

        HighPrecision integralResult = IntegrationUtils::integrate(integrand, this->p.t0, t);
        HighPrecision result = prefactor * integralResult;

        {
            std::lock_guard<std::mutex> lock(cacheMutex);
            rhoPhiCache[key] = result;
        }

        return result;
    };
}


EnergyDensity PhiParticle::energyDensityMatter(HighPrecision t0)
{
    constexpr double n = 4.0;
    ChiDecayRate chiDecay(this->p, n, t0);
    return [this, t0, chiDecay](HighPrecision t)->HighPrecision{
        // n = 4 in matter dominated Universe
        return pow((t0 / t), HighPrecision(2.0)) * exp(-chiDecay(t)) * this->getInitialRhoMatter();
    };
}


EnergyDensity PhiParticle::energyDensityRadiation(HighPrecision t0)
{
    constexpr double n = 2.0;
    ChiDecayRate chiDecay(this->p, n, t0);
    return [this, t0, chiDecay](HighPrecision t)->HighPrecision{
        // n = 2 in matter dominated Universe
        return pow((t0 / t), HighPrecision(3.0 / 2.0)) * exp(-chiDecay(t)) * this->getInitialRhoRadiation();
    };    
}