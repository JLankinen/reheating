#ifndef PHI_PARTICLE_H_
#define PHI_PARTICLE_H_

#include <map>
#include <mutex>
#include <optional>

#include "parameters/parameters.hpp"
#include "utils/types.hpp"

/**
 * @brief Represents the massive Phi particle in the model.
 * 
 * This class traces the evolution of a massive Phi particle and provides methods related to the
 * energy density evolution of the massive phi particle. 
 * 
 * Returns an energy density function which can be used in solvers
 * or further analysis.
 */
class PhiParticle
{
    private:
        ModelParameters p;
        std::map<double, double> rhoPhiCache;
        std::mutex cacheMutex;
        double initialRhoMatter;
        double initialRhoRadiation;
    public:
        explicit PhiParticle(const ModelParameters& _p) : p{_p} {};
        double creationRate(double t);
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