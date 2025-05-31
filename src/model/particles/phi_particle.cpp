#include <cmath>

#include "model/particles/phi_particle.hpp"
#include "model/energy/creation_decay.hpp"
#include "utils/types.hpp"
#include "utils/integration.hpp"


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

HighPrecision PhiParticle::creationRate(HighPrecision t)
{
    HighPrecision arg = -pow((HighPrecision(3.0) * p.m * t / HighPrecision(2.0)), HighPrecision(2.0) / HighPrecision(3.0));

    HighPrecision airyAi = boost::math::airy_ai(arg);
    HighPrecision airyBi = boost::math::airy_bi(arg);

    HighPrecision airySum = pow(airyAi, 2.0) + pow(airyBi, 2.0);

    HighPrecision factor = HighPrecision(3.0) * pow(p.m * p.b, HighPrecision(13.0) / HighPrecision(3.0)) / (HighPrecision(32.0) * p.b);

    return factor * t * airySum;
}

EnergyDensity PhiParticle::energyDensityStiff()
{
    constexpr double n = 1.0;
    ChiDecayRate chiDecay(this->p, n, this->p.t0);

    return [this, chiDecay](HighPrecision t) -> HighPrecision
    {

       HighPrecision key = quantize(t, 30);  // Quantize to 30 digits for cache key

        {
            std::lock_guard<std::mutex> lock(this->cacheMutex);
            auto it = this->rhoPhiCache.find(key);
            if (it != this->rhoPhiCache.end()) {
                return it->second;
            }
        }
        
        HighPrecision prefactor = (HighPrecision(1.0) / t) * exp(-chiDecay(t));

        auto integrand = [&](HighPrecision tprime)
        {
            HighPrecision val = tprime * this->creationRate(tprime) * exp(chiDecay(tprime));
            return val;
        };

        HighPrecision integralResult = IntegrationUtils::integrate(integrand, this->p.t0, t);
        HighPrecision result = prefactor * integralResult;

        {
            std::lock_guard<std::mutex> lock(this->cacheMutex);
            this->rhoPhiCache[key] = result;
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