#include <cmath>

#include "model/particles/chi_particle.hpp"
#include "model/particles/phi_particle.hpp"
#include "model/energy/creation_decay.hpp"
#include "utils/integration.hpp"


void ChiParticle::setInitialRhoMatter(const HighPrecision& rhoInit)
{
    this->initialRhoMatter = rhoInit;
}

HighPrecision ChiParticle::getInitialRhoMatter() const
{
    return initialRhoMatter;
}

void ChiParticle::setInitialRhoRadiation(const HighPrecision& rhoInit)
{
    this->initialRhoRadiation = rhoInit;
}

HighPrecision ChiParticle::getInitialRhoRadiation() const
{
    return initialRhoRadiation;
}

EnergyDensity ChiParticle::energyDensityStiff()
{   
    constexpr double n = 1.0; // For stiff matter universe, n = 1.
    ChiDecayRate chiDecay(this->p, n, this->p.t0);
    return [this, chiDecay](HighPrecision t) -> HighPrecision
    {
        HighPrecision prefactor = 1 / pow(t, HighPrecision(4.0 / 3.0));
        EnergyDensity rhoPhi = phiParticle->energyDensityStiff();
        
        auto integrand = [&](HighPrecision tprime)
        {
            HighPrecision val = chiDecay(tprime) * rhoPhi(tprime)
                                * pow(tprime, HighPrecision(4.0) / HighPrecision(3.0));
            return val;
        };
    
        HighPrecision integralResult = IntegrationUtils::integrate(integrand, this->p.t0, t);
        return prefactor * integralResult;
    };
}


EnergyDensity ChiParticle::energyDensityMatter(HighPrecision t0)
{
    constexpr double n = 4.0; // n=4 for matter domination 
    ChiDecayRate chiDecay(this->p, n, t0);

    return [this, t0, chiDecay](HighPrecision t)->HighPrecision{
        HighPrecision initialRho = this->getInitialRhoMatter() * pow(t0 / t, HighPrecision(8.0 / 3.0));
        HighPrecision prefactor = pow(t, - HighPrecision(8.0 / 3.0));
        EnergyDensity rhoMat = phiParticle->energyDensityMatter(t0);
        auto integrand = [&] (HighPrecision tprime)
        {
            HighPrecision val = chiDecay(tprime) * rhoMat(tprime) * pow(tprime, (8.0 / 3.0));
            return val;
        };

        auto integral = IntegrationUtils::integrate(integrand, t0, t);
        return prefactor * integral + initialRho;
    };
}

EnergyDensity ChiParticle::energyDensityRadiation(HighPrecision t0)
{
    constexpr double n = 2.0; // n = 2 in radiation dominated universe
    ChiDecayRate chiDecay(this->p, n, t0);

    return [this, t0, chiDecay](HighPrecision t)->HighPrecision{
        HighPrecision time = pow(t0, (1.0 / 2.0)) / pow(t, (1.0 / 2.0));
        HighPrecision initialRho = this->getInitialRhoMatter() * pow(time, 4);  // Rho_chi_mat(tau_eq)
        HighPrecision prefactor = 1 / pow(t, 2);
        EnergyDensity rhoRad = phiParticle->energyDensityRadiation(t0);

        auto integrand = [&] (HighPrecision tprime)
        {
            HighPrecision val = chiDecay(tprime) * rhoRad(tprime) * pow(tprime, 2.0);
            return val;
        };

        auto integral = IntegrationUtils::integrate(integrand, t0, t);
        return prefactor * integral + initialRho;
    };
}

