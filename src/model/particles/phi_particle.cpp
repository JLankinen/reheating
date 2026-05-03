#include <cmath>
#include <sstream>
#include <iomanip>

#include "model/particles/phi_particle.hpp"
#include "model/energy/creation_decay.hpp"
#include "utils/types.hpp"
#include "utils/integration.hpp"


double quantize(const double& value, int digits = 30) {
    std::ostringstream oss;
    oss << std::setprecision(digits) << std::fixed << value;
    return static_cast<double>(std::stod(oss.str()));
}

void PhiParticle::setInitialRhoMatter(const double& rhoInit)
{
    this->initialRhoMatter = rhoInit;
}

double PhiParticle::getInitialRhoMatter() const
{
    return initialRhoMatter;
}

void PhiParticle::setInitialRhoRadiation(const double& rhoInit)
{
    this->initialRhoRadiation = rhoInit;
}

double PhiParticle::getInitialRhoRadiation() const
{
    return initialRhoRadiation;
}

double PhiParticle::creationRate(double t)
{
    double arg = -pow((3.0 * p.m * t / 2.0), 2.0 / 3.0);

    double airyAi = boost::math::airy_ai(arg);
    double airyBi = boost::math::airy_bi(arg);

    double airySum = pow(airyAi, 2.0) + pow(airyBi, 2.0);

    double factor = 3.0 * pow(p.m * p.b, 13.0 / 3.0) / (32.0 * p.b);

    return factor * t * airySum;
}

EnergyDensity PhiParticle::energyDensityStiff()
{
    constexpr double n = 1.0;
    //constexpr double n = 0.0;
    ChiDecayRate chiDecay(this->p, n, this->p.t0);

    return [this, chiDecay](double t) -> double
    {

       double key = quantize(t, 30);  // Quantize to 30 digits for cache key

        {
            std::lock_guard<std::mutex> lock(this->cacheMutex);
            auto it = this->rhoPhiCache.find(key);
            if (it != this->rhoPhiCache.end()) {
                return it->second;
            }
        }
        
        double prefactor = (1.0 / t) * exp(-chiDecay(t));

        auto integrand = [&](double tprime)
        {
            double val = tprime * this->creationRate(tprime) * exp(chiDecay(tprime));
            return val;
        };

        double integralResult = IntegrationUtils::integrate(integrand, this->p.t0, t);
        double result = prefactor * integralResult;

        {
            std::lock_guard<std::mutex> lock(this->cacheMutex);
            this->rhoPhiCache[key] = result;
        }


        return result;
    };
}


EnergyDensity PhiParticle::energyDensityMatter(double t0)
{
    constexpr double n = 4.0;
    //constexpr double n = 0.0;
    ChiDecayRate chiDecay(this->p, n, t0);
    return [this, t0, chiDecay](double t)->double{
        // n = 4 in matter dominated Universe
        return pow((t0 / t), 2.0) * exp(-chiDecay(t)) * this->getInitialRhoMatter();
    };
}


EnergyDensity PhiParticle::energyDensityRadiation(double t0)
{
    constexpr double n = 2.0;
    //constexpr double n = 0.0;
    ChiDecayRate chiDecay(this->p, n, t0);
    return [this, t0, chiDecay](double t)->double{
        // n = 2 in radiation dominated Universe
        return pow((t0 / t), 3.0 / 2.0) * exp(-chiDecay(t)) * this->getInitialRhoRadiation();
    };    
}