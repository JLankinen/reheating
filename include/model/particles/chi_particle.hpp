#ifndef CHI_PARTICLE_H_
#define CHI_PARTICLE_H_

#include "parameters/parameters.hpp"
#include "utils/types.hpp"

/**
 * @brief Represents the massless chi particle in the model.
 * 
 * This class provides methods related to the energy density evolution
 * of the massless phi particle (radiation). 
 * 
 * Returns an energy density function which can be used in solvers
 * or further analysis.
 */
class ChiParticle
{
    private:
        ModelParameters p;
        EnergyDensity rhoPhi;
    public:
        explicit ChiParticle(const ModelParameters& _p, EnergyDensity _rhoPhi) :
        p{_p}, rhoPhi{_rhoPhi} {};
  
        // Returns RhoChiStiff as a function of t.
        EnergyDensity energyDensityStiff();
};


#endif