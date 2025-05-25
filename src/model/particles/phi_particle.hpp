#ifndef PHI_PARTICLE_H_
#define PHI_PARTICLE_H_

#include "parameters/parameters.hpp"
#include "utils/types.hpp"

class PhiParticle
{
    private:
        ModelParameters p;
    
    public:
        explicit PhiParticle(const ModelParameters& _p) : p{_p} {};
        
        // Returns RhoPhiStiff as a function of t.
        EnergyDensity energyDensityStiff();
};

#endif