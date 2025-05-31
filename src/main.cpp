#include <iostream>
#include <iomanip>
#include <chrono>
#include <memory>

#include "parameters/parameters.hpp"
#include "solvers/equal_time_solver.hpp"
#include "model/particles/phi_particle.hpp"
#include "model/particles/chi_particle.hpp"
#include "model/particles/stiff_matter.hpp"
#include "simulation/simulation.hpp"
#include "utils/maximization.hpp"
#include "utils/integration.hpp"

int main()
{
    /**TODO
     * 1. Fix the restriction
     * (2. Add refinement e.g., newton to solver builder)
     * 3. Fix numerical issues
     *   3.1 Find a way to bracket the upper automatically tight.
     *    Make custom bracketer: increase by 1e1 until sign change.
     * 4. Add logging
     * 5. Create simulation object. Multithread to each core.
     * 6. Add output to parquet
     * 7. Generate parameter space.
     */
    ModelParameters p;
    p.logParameters();

    

    // Make one shared PhiParticle to be used throughout the program.
    auto phiParticle = std::make_shared<PhiParticle>(p);
    auto chiParticle = ChiParticle(p, phiParticle);
    auto stiff = StiffMatter(p);

    // Find equal time starting from the stiff matter phase.
    EnergyDensity rhoPhiStiff = phiParticle->energyDensityStiff();
    EnergyDensity rhoStiff = stiff.energyDensity();
    EnergyDensity rhoChiStiff = chiParticle.energyDensityStiff();

    auto stiffSolver = EqualTimeSolver(rhoStiff, rhoPhiStiff, p.t0).withRestriction(rhoChiStiff);
    stiffSolver.setUpperLimit(p.t0 * HighPrecision("1e27"));
    auto [t_eq, rhoStiffEq, rhoPhiEq] = stiffSolver.getEqualTime();

    std::cout << "t_eq: " << t_eq << " RhoPhi at t_eq: " << rhoPhiEq << " RhoStiff at t_eq: " << rhoStiffEq << std::endl;
        
    // Use t_eq as initial time for the calculations in the matter dominated phase.
    // Set the initial values at t_eq
    HighPrecision rhoChiEq = rhoChiStiff(t_eq);
    phiParticle->setInitialRhoMatter(rhoPhiEq);
    chiParticle.setInitialRhoMatter(rhoChiEq);
    std::cout << "Initial RhoChi end of stiff: " << rhoChiEq << std::endl;
    std::cout << "Initial RhoPhi end of stiff: " << rhoPhiEq << std::endl;

    EnergyDensity rhoPhiMat = phiParticle->energyDensityMatter(t_eq);
    EnergyDensity rhoChiMat = chiParticle.energyDensityMatter(t_eq);
    auto radMatSolver = EqualTimeSolver(rhoPhiMat, rhoChiMat, t_eq);
    radMatSolver.setUpperLimit(t_eq * HighPrecision("1e18"));

    auto [tau_eq, rhoPhiMatEq, rhoChiMatEq] = radMatSolver.getEqualTime();

    std::cout << "tau_eq: " << tau_eq << " RhoPhi at tau_eq: " << rhoPhiMatEq << " RhoChi at tau_eq: " << rhoChiMatEq << std::endl;

    // Use tau_eq as initial time for the calculations in the radiation dominated phase.
    phiParticle->setInitialRhoRadiation(rhoPhiMatEq);
    chiParticle.setInitialRhoRadiation(rhoChiMatEq);
    std::cout << "Initial RhoChi end of matter: " << rhoPhiMatEq << std::endl;
    std::cout << "Initial RhoPhi end of matter: " << rhoChiMatEq << std::endl;

    EnergyDensity rhoPhiRad = phiParticle->energyDensityRadiation(tau_eq);
    EnergyDensity rhoChiRad = chiParticle.energyDensityRadiation(tau_eq);
    auto radSolver = EqualTimeSolver(rhoPhiRad, rhoChiRad, tau_eq);
    radSolver.setUpperLimit(tau_eq * HighPrecision("1e10"));
    auto [tau2_eq, rhoPhiRadEq, rhoChiRadEq] = radSolver.getEqualTime();

    std::cout << "tau2_eq: " << tau2_eq << " RhoPhi at tau2_eq: " << rhoPhiRadEq << " RhoChi at tau2_eq: " << rhoChiRadEq << std::endl;

    // Maximize rhochirad
    auto t_rh = maximize(rhoChiRad, tau2_eq, tau2_eq * HighPrecision("1e5"));
    HighPrecision reheatingTemperature = IntegrationUtils::integrate(chiParticle.energyDensityRadiation(tau2_eq), tau2_eq, t_rh);
    std::cout << "Reheating time: " << t_rh << ", Reheating temperature: " << pow(reheatingTemperature, HighPrecision(1.0 / 4.0)) << std::endl;
    
}