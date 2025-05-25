#ifndef CHI_PARTICLE_H_
#define CHI_PARTICLE_H_

#include "parameters/parameters.hpp"
#include "utils/types.hpp"

class ChiParticle
{
    private:
        ModelParameters p;
        EnergyDensity rhoPhi;
    public:
        explicit ChiParticle(const ModelParameters& _p, EnergyDensity _rhoPhi) :
        p{_p}, rhoPhi{_rhoPhi} {};
  
        // Returns RhoPhiStiff as a function of t.
        EnergyDensity energyDensityStiff();
};


#endif