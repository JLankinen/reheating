#ifndef CHI_PARTICLE_H_
#define CHI_PARTICLE_H_

#include <memory>

#include "parameters/parameters.hpp"
#include "utils/types.hpp"
#include "model/particles/phi_particle.hpp"

/**
 * @brief Represents the massless chi particle in the model.
 * 
 * This class traces the evolution of a massless Chi particle and provides methods related to the
 * energy density evolution of the massless chi particle. 
 * 
 * As a decay product of the massive Phi particle, this particle class is coupled to phi class.
 * 
 * Returns an energy density function which can be used in solvers
 * or further analysis.
 */
class ChiParticle
{
    private:
        ModelParameters p;
        std::shared_ptr<PhiParticle> phiParticle;
        double initialRhoMatter;
        double initialRhoRadiation;
    public:
        explicit ChiParticle(const ModelParameters& _p, std::shared_ptr<PhiParticle> phi) :
        p{_p}, phiParticle(std::move(phi)) {};
  
        // Returns RhoChiStiff as a function of t.
        EnergyDensity energyDensityStiff();
        EnergyDensity energyDensityMatter(double t0);
        EnergyDensity energyDensityRadiation(double t0);
        // Getters and Setters to set initial values once obtained
        void setInitialRhoMatter(const double& rhoInit);
        double getInitialRhoMatter() const;
        void setInitialRhoRadiation(const double& rhoInit);
        double getInitialRhoRadiation() const;
};


#endif