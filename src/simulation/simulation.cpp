#include "model/energy/creation_decay.hpp"
#include "model/particles/chi_particle.hpp"
#include "model/particles/phi_particle.hpp"
#include "model/particles/stiff_matter.hpp"
#include "solvers/equal_time_solver.hpp"
#include "simulation/simulation.hpp"
#include "parameters/parameters.hpp"
#include "utils/maximization.hpp"
#include "utils/integration.hpp"


Simulation::Simulation(const ModelParameters& p_) :
    p{p_},
    phi{std::make_shared<PhiParticle>(p_)},
    chi{ChiParticle(p_, phi)},
    stiff{StiffMatter(p_)}
    {};


SimulationResults Simulation::run()
{
    // Return time of equality and energy densities of stiff matter and that particle
    // (massive phi or massles chi) depending on which one reaches equality first.
    auto [toMatter, t_eq, rhoStiffEq, rhoEq, bothFound] = runStiffPhase();

    if (toMatter)
    {
        // Set the initial values
        EnergyDensity rhoChiStiff = chi.energyDensityStiff();
        HighPrecision rhoChiEq = rhoChiStiff(t_eq);
        phi->setInitialRhoMatter(rhoEq);
        chi.setInitialRhoMatter(rhoChiEq);
        
        auto [tau_eq, rhoPhiMatEq, rhoChiMatEq] = runMatterPhase(t_eq);
        phi->setInitialRhoRadiation(rhoPhiMatEq);
        chi.setInitialRhoRadiation(rhoChiMatEq);
        
        //auto [t_eq_rad, rhoPhiRadEq, rhoChiRadEq] = runRadiationPhase(tau_eq);
        auto [tempRH, timeRH] = getReheatingTemperatureAndTime(tau_eq);
    
        return SimulationResults{
            .params = p,
            .reheating_temp = tempRH,
            .reheating_time = timeRH,
            .t_eq = t_eq,
            .rhoStiff_t_eq = rhoStiffEq,
            .rhoPhiStiff_t_eq = rhoEq,
            .rhoChi_t_eq = rhoChiEq,
            .tau_eq = tau_eq,
            .rhoPhiMatEq = rhoPhiMatEq,
            .rhoChiMatEq = rhoChiMatEq,
            .toMatter = true,
            .bothFound = bothFound 
            };
    }
    else
    {
        // Now rhoEq is the value of massless particles (rho_chi) at t_eq.
        // Calculate value of massive particles at t_eq.
        EnergyDensity rhoPhiStiff = phi->energyDensityStiff();
        HighPrecision rhoPhiEq = rhoPhiStiff(t_eq);

        // Set initial conditions
        phi->setInitialRhoRadiation(rhoPhiEq);
        chi.setInitialRhoRadiation(rhoEq);
        auto [tempRH, timeRH] = getReheatingTemperatureAndTime(t_eq);

        return SimulationResults{
        .params = p,
        .reheating_temp = tempRH,
        .reheating_time = timeRH,
        .t_eq = t_eq,
        .rhoStiff_t_eq = rhoStiffEq,
        .rhoChi_t_eq = rhoEq,
        .tau_eq = 0,
        .rhoPhiMatEq = 0,
        .rhoChiMatEq = 0,
        .toMatter = false,
        .bothFound = bothFound
        };
    }

}


/**
 * @brief Check if the energy density of massive particles is greater than massless particles
 * at the time of equality. If it is, Universe ends up in matter domination. Otherwise it ends
 * up in radiation domination.
 */
bool Simulation::toMatter(const EnergyDensity &rhoChi, HighPrecision rhoPhi, HighPrecision timeEquality)
{
    return rhoPhi > rhoChi(timeEquality);
}


std::tuple<bool, HighPrecision, HighPrecision, HighPrecision, bool> Simulation::runStiffPhase()
{
    EnergyDensity rhoChiStiff = chi.energyDensityStiff();
    EnergyDensity rhoPhiStiff = phi->energyDensityStiff();
    EnergyDensity rhoStiff = stiff.energyDensity();

    // Compute the equality times if they exist
    auto stiffPhi = EqualTimeSolver(rhoStiff, rhoPhiStiff, p.t0).findEqualTime();  // Equal time for stiff and phi
    auto stiffChi = EqualTimeSolver(rhoStiff, rhoChiStiff, p.t0).findEqualTime();  // Equal time for stiff and chi

    // Pick the one which is earlier
    if (stiffPhi && stiffChi)
    {
        auto [t_phi, rhoStiffPhi, RhoPhiEq] = *stiffPhi;
        auto [t_chi, rhoStiffChi, RhoChiEq] = *stiffChi;

        if (t_phi < t_chi)
        {
            //std::cout << "Both found -> to matter" << std::endl;
            // First bool: toMatter, last bool: both found.
            return {true, t_phi, rhoStiffPhi, RhoPhiEq, true};  // -> matter
        }
        else
        {
            //std::cout << "Both found -> to radiation" << std::endl;
            return {false, t_chi, rhoStiffChi, RhoChiEq, true}; // -> radiation
        }
    }

    if (stiffPhi)
    {
        auto [t_phi, rhoStiffEq, RhoPhiEq] = *stiffPhi;
        //std::cout << "Only phi found --> matter" << std::endl;
        return {true, t_phi, rhoStiffEq, RhoPhiEq, false};  // -> matter
    }

    if (stiffChi)
    {
        auto [t_chi, rhoStiffEq, RhoChiEq] = *stiffChi;
        //std::cout << "Only chi found --> radiation" << std::endl;
        return {false, t_chi, rhoStiffEq, RhoChiEq, false}; // -> radiation
    }

    // Something is wrong if this is hit.
    throw std::runtime_error("No transition from stiff phase.");
}



std::tuple<HighPrecision, HighPrecision, HighPrecision> Simulation::runMatterPhase(HighPrecision t0)
{
    EnergyDensity rhoPhiMat = phi->energyDensityMatter(t0);
    EnergyDensity rhoChiMat = chi.energyDensityMatter(t0);
    auto radMatSolver = EqualTimeSolver(rhoPhiMat, rhoChiMat, t0);

    auto [tau_eq, rhoPhiMatEq, rhoChiMatEq] = radMatSolver.getEqualTime();

    return {tau_eq, rhoPhiMatEq, rhoChiMatEq};
}

std::tuple<HighPrecision, HighPrecision, HighPrecision> Simulation::runRadiationPhase(HighPrecision t0)
{
    EnergyDensity rhoPhiRad = phi->energyDensityRadiation(t0);
    EnergyDensity rhoChiRad = chi.energyDensityRadiation(t0);
    auto radSolver = EqualTimeSolver(rhoPhiRad, rhoChiRad, t0);
    auto [t_eq_rad, rhoPhiRadEq, rhoChiRadEq] = radSolver.getEqualTime();

    return {t_eq_rad, rhoPhiRadEq, rhoChiRadEq};
}

std::pair<HighPrecision, HighPrecision> Simulation::getReheatingTemperatureAndTime(HighPrecision tau_eq)
{
    EnergyDensity rhoChiRad = this->chi.energyDensityRadiation(tau_eq);
    auto t_rh = maximize(rhoChiRad, tau_eq, tau_eq * HighPrecision("1e5"));
    HighPrecision reheatingTemperature = IntegrationUtils::integrate(this->chi.energyDensityRadiation(tau_eq),tau_eq, t_rh);
    HighPrecision T_RH = pow(reheatingTemperature, HighPrecision(1.0 / 4.0));
    return std::pair(T_RH, t_rh);
}
