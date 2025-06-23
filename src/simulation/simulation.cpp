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
    p.logParameters();
    auto [t_eq, rhoStiffEq, rhoPhiEq] = runStiffPhase();
    auto [tau_eq, rhoPhiMatEq, rhoChiMatEq] = runMatterPhase(t_eq);
    auto [t_eq_rad, rhoPhiRadEq, rhoChiRadEq] = runRadiationPhase(tau_eq);
    auto [tempRH, timeRH] = getReheatingTemperatureAndTime(tau_eq, t_eq_rad);

    return SimulationResults{
        .params = p,
        .reheating_temp = tempRH,
        .reheating_time = timeRH,
        .t_eq = t_eq,
        .rhoStiff_t_eq = rhoStiffEq,
        .rhoPhiStiff_t_eq = rhoPhiEq,
        .tau_eq = tau_eq,
        .rhoPhiMatEq = rhoPhiMatEq,
        .rhoChiMatEq = rhoChiMatEq,
        .tau2_eq = t_eq_rad,
        .rhoPhiRadEq = rhoPhiRadEq,
        .rhoChiRadEq = rhoChiRadEq  
    };
}


std::tuple<HighPrecision, HighPrecision, HighPrecision> Simulation::runStiffPhase()
{
        /* Find equal time starting from the stiff matter phase. */
        EnergyDensity rhoPhiStiff = phi->energyDensityStiff();
        EnergyDensity rhoStiff = stiff.energyDensity();
        EnergyDensity rhoChiStiff = chi.energyDensityStiff();
    
        auto stiffSolver = EqualTimeSolver(rhoStiff, rhoPhiStiff, p.t0).withRestriction(rhoChiStiff);
        auto [t_eq, rhoStiffEq, rhoPhiEq] = stiffSolver.getEqualTime();
    
        std::cout << "t_eq: " << t_eq << " RhoPhi at t_eq: " << rhoPhiEq << " RhoStiff at t_eq: " << rhoStiffEq << std::endl;

        // Set the initial values
        HighPrecision rhoChiEq = rhoChiStiff(t_eq);
        phi->setInitialRhoMatter(rhoPhiEq);
        chi.setInitialRhoMatter(rhoChiEq);
        std::cout << "Initial RhoChi end of stiff: " << rhoChiEq << std::endl;
        std::cout << "Initial RhoPhi end of stiff: " << rhoPhiEq << std::endl;
        return std::make_tuple(t_eq, rhoStiffEq, rhoPhiEq);
}

std::tuple<HighPrecision, HighPrecision, HighPrecision> Simulation::runMatterPhase(HighPrecision t0)
{
    EnergyDensity rhoPhiMat = phi->energyDensityMatter(t0);
    EnergyDensity rhoChiMat = chi.energyDensityMatter(t0);
    auto radMatSolver = EqualTimeSolver(rhoPhiMat, rhoChiMat, t0);

    auto [tau_eq, rhoPhiMatEq, rhoChiMatEq] = radMatSolver.getEqualTime();

    std::cout << "tau_eq: " << tau_eq << " RhoPhi at tau_eq: " << rhoPhiMatEq << " RhoChi at tau_eq: " << rhoChiMatEq << std::endl;

    /* Use tau_eq as initial time for the calculations in the radiation dominated phase.*/
    phi->setInitialRhoRadiation(rhoPhiMatEq);
    chi.setInitialRhoRadiation(rhoChiMatEq);
    std::cout << "Initial RhoChi end of matter: " << rhoPhiMatEq << std::endl;
    std::cout << "Initial RhoPhi end of matter: " << rhoChiMatEq << std::endl;

    return std::make_tuple(tau_eq, rhoPhiMatEq, rhoChiMatEq);
}

std::tuple<HighPrecision, HighPrecision, HighPrecision> Simulation::runRadiationPhase(HighPrecision t0)
{
    EnergyDensity rhoPhiRad = phi->energyDensityRadiation(t0);
    EnergyDensity rhoChiRad = chi.energyDensityRadiation(t0);
    auto radSolver = EqualTimeSolver(rhoPhiRad, rhoChiRad, t0);
    auto [t_eq_rad, rhoPhiRadEq, rhoChiRadEq] = radSolver.getEqualTime();

    std::cout << "t_eq rad: " << t_eq_rad << " RhoPhi at t_eq: " << rhoPhiRadEq << " RhoChi at t_eq: " << rhoChiRadEq << std::endl;

    return std::make_tuple(t_eq_rad, rhoPhiRadEq, rhoChiRadEq);
}

std::pair<HighPrecision, HighPrecision> Simulation::getReheatingTemperatureAndTime(HighPrecision tau_eq, HighPrecision tau2_eq)
{
    EnergyDensity rhoChiRad = this->chi.energyDensityRadiation(tau_eq);
    auto t_rh = maximize(rhoChiRad, tau2_eq, tau2_eq * HighPrecision("1e5"));
    HighPrecision reheatingTemperature = IntegrationUtils::integrate(this->chi.energyDensityRadiation(tau2_eq), tau2_eq, t_rh);
    std::cout << "Reheating time: " << t_rh << ", Reheating temperature: " << pow(reheatingTemperature, HighPrecision(1.0 / 4.0)) << std::endl;

    return std::pair(reheatingTemperature, t_rh);
}
