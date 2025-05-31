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
        std::map<HighPrecision, HighPrecision> rhoPhiCache;
        std::mutex cacheMutex;
        HighPrecision initialRhoMatter;
        HighPrecision initialRhoRadiation;
    public:
        explicit PhiParticle(const ModelParameters& _p) : p{_p} {};
        HighPrecision creationRate(HighPrecision t);
        EnergyDensity energyDensityStiff();
        EnergyDensity energyDensityMatter(HighPrecision t0);
        EnergyDensity energyDensityRadiation(HighPrecision t0);
        // Getters and Setters to set initial values once obtained
        void setInitialRhoMatter(const HighPrecision& rhoInit);
        HighPrecision getInitialRhoMatter() const;
        void setInitialRhoRadiation(const HighPrecision& rhoInit);
        HighPrecision getInitialRhoRadiation() const;



};

#endif