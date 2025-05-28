#ifndef PHI_PARTICLE_H_
#define PHI_PARTICLE_H_

#include <map>
#include <mutex>
#include "parameters/parameters.hpp"
#include "utils/types.hpp"

/**
 * @brief Represents the massive Phi particle in the model.
 * 
 * This class provides methods related to the energy density evolution
 * of the massive phi particle. 
 * 
 * Returns an energy density function which can be used in solvers
 * or further analysis.
 */
class PhiParticle
{
    private:
        ModelParameters p;
        static std::map<HighPrecision, HighPrecision> rhoPhiCache;
        static std::mutex cacheMutex;
    public:
        explicit PhiParticle(const ModelParameters& _p) : p{_p} {};

        EnergyDensity energyDensityStiff();
        EnergyDensity energyDensityMatter(HighPrecision t0, HighPrecision intialRho);
};

#endif