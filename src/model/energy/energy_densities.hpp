#ifndef ENERGY_DENSITY_H_
#define ENERGY_DENSITY_H_

#include "parameters/parameters.hpp"
#include "model/energy/creation_decay.hpp"

/**
 * @brief Calculates the stiff matter energy density rho_stiff at time t$.
 * 
 * @param p Model parameters containing gravitational constant and initial conditions
 * @param t Time variable at which to evaluate the energy density
 * @return HighPrecision Energy density rho_stiff(t)
 */
HighPrecision RhoStiff(ModelParameters p, HighPrecision t);

/**
 * @brief Computes the energy density rho_phi related to the massive scalar field phi decaying
 * into massless chi particles in stiff matter dominated universe.
 * 
 * Evaluates rho_phi at time t by numerically integrating a decay-influenced function,
 * including caching for performance.
 * 
 * Uses Gauss-Kronrod quadrature for numerical integration within the interval [p.t0, t].
 * 
 * @param p Model parameters defining decay rates, initial time, and other physics constants
 * @param t Time variable at which to evaluate the energy density
 * @return HighPrecision Energy density rho_phi(t)
 */
HighPrecision RhoPhiStiff(ModelParameters p, HighPrecision t);

/**
 * @brief Calculates the energy density rho_chi of the massless chi particles.
 * 
 * Computes rho_chi(t) through numerical integration involving decay rates and rho_phi.
 * 
 * The function integrates over the interval [p.t0, t] using Gauss-Kronrod quadrature.
 * 
 * @param p Model parameters defining the decay and initial conditions
 * @param t Time variable at which to evaluate the energy density
 * @return HighPrecision Energy density rho_chi(t)
 */
HighPrecision RhoChiStiff(ModelParameters p, HighPrecision t);

/**
 * @brief Finds the time at which two energy densities are equal, subject to a constraint.
 *
 * Solves for the time `t` such that `f1(p, t) == f2(p, t)` in the interval [`lowerLimit`, `upperLimit`],
 * under the constraint that `f2(p, t) > restrictionFunc(p, t)`.
 *
 * @param f1 First function (e.g., energy density of phi)
 * @param f2 Second function (e.g., energy density of chi)
 * @param restrictionFunc Constraint function (e.g., radiation or stiff background bound)
 * @param p Model parameters
 * @param lowerLimit Lower bound of the search interval
 * @param upperLimit Upper bound of the search interval
 * @return A pair containing:
 *         - The time `t_eq` where `f1(p, t) == f2(p, t)`
 *         - The value of `f(p, t_eq)` at that time
 * @throws std::runtime_error if no root satisfies the constraint
 */
std::pair<HighPrecision, HighPrecision> EqualTime(std::function<HighPrecision(ModelParameters, HighPrecision t)> &rhoStiff,
                 std::function<HighPrecision(ModelParameters, HighPrecision t)> &rhoPhiStiff,
                 std::function<HighPrecision(ModelParameters, HighPrecision t)> restrictionFunc,
                 ModelParameters p,
                 HighPrecision lowerLimit,
                 HighPrecision upperLimit);

#endif