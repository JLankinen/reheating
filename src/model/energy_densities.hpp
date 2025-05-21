#ifndef ENERGY_DENSITY_H_
#define ENERGY_DENSITY_H_

#include "parameters/parameters.hpp"
#include "model/creation_decay.hpp"

/**
 * Stiff matter energy density
 */
HighPrecision RhoStiff(HighPrecision t);

/**
 * Energy density for massive phi particles during stiff matter era.
 */
HighPrecision RhoPhiStiff(ModelParameters p, HighPrecision t);

/**
 * Energy density for massless chi particles during stiff matter era.
 */
HighPrecision RhoChiStiff(ModelParameters p, HighPrecision t);

/**
 * Calculate t_eq for stiff matter and massive particle energy densities.
 */
HighPrecision EqualTime(std::function<HighPrecision(HighPrecision t)> rhoStiff,
                 std::function<HighPrecision(ModelParameters, HighPrecision t)> rhoPhiStiff,
                 std::function<HighPrecision(ModelParameters, HighPrecision t)> restrictionFunc,
                 ModelParameters p,
                 HighPrecision lowerLimit,
                 HighPrecision upperLimit);

#endif