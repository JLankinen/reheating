#ifndef STIFF_MATTER_H_
#define STIFF_MATTER_H_

#include "parameters/parameters.hpp"
#include "utils/types.hpp"

/**
 * @brief Represents a stiff matter background in the cosmological model.
 *
 * Provides a callable function that returns the energy density 
 * as a function of time.
 */
class StiffMatter final
{
    private:
        ModelParameters p;
    public:
        explicit StiffMatter(const ModelParameters& _p) : p{_p} {};
    /**
     * @brief Returns a callable function representing the energy density over time.
     * 
     * The energy density is computed as:
     *     ρ(t) = 1 / (24 * π * G_N * t²)
     * 
     * @return EnergyDensity function (HighPrecision -> HighPrecision)
     */
        EnergyDensity energyDensity();
};

#endif