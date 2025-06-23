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

/**
 * @brief Check if the energy density of massive particles is greater than massless particles
 * at the time of equality. If it is, Universe ends up in matter domination. Otherwise it ends
 * up in radiation domination.
 */
bool Simulation::toMatter(HighPrecision rhoPhi, HighPrecision timeEquality)
{
    EnergyDensity rhoChiStiff = chi.energyDensityStiff();
    return rhoPhi > rhoChiStiff(timeEquality);
}


SimulationResults Simulation::run()
{
    p.logParameters();
    auto [t_eq, rhoStiffEq, rhoPhiEq] = runStiffPhase();

    if (toMatter(rhoPhiEq, t_eq))
    {
        // Set the initial values
        EnergyDensity rhoChiStiff = chi.energyDensityStiff();
        HighPrecision rhoChiEq = rhoChiStiff(t_eq);
        phi->setInitialRhoMatter(rhoPhiEq);
        chi.setInitialRhoMatter(rhoChiEq);
        
        auto [tau_eq, rhoPhiMatEq, rhoChiMatEq] = runMatterPhase(t_eq);
        phi->setInitialRhoRadiation(rhoPhiMatEq);
        chi.setInitialRhoRadiation(rhoChiMatEq);
        
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
            .rhoChiRadEq = rhoChiRadEq,
            .toMatter = true              
            };
    }
    else
    {
        EnergyDensity rhoChiStiff = chi.energyDensityStiff();
        HighPrecision rhoChiEq = rhoChiStiff(t_eq);
        phi->setInitialRhoRadiation(rhoPhiEq);
        chi.setInitialRhoRadiation(rhoChiEq);
        auto [t_eq_rad, rhoPhiRadEq, rhoChiRadEq] = runRadiationPhase(t_eq);
        auto [tempRH, timeRH] = getReheatingTemperatureAndTime(t_eq, t_eq_rad);

        return SimulationResults{
        .params = p,
        .reheating_temp = tempRH,
        .reheating_time = timeRH,
        .t_eq = t_eq,
        .rhoStiff_t_eq = rhoStiffEq,
        .rhoPhiStiff_t_eq = rhoPhiEq,
        .tau_eq = 0,
        .rhoPhiMatEq = 0,
        .rhoChiMatEq = 0,
        .tau2_eq = t_eq_rad,
        .rhoPhiRadEq = rhoPhiRadEq,
        .rhoChiRadEq = rhoChiRadEq,
        .toMatter = false
        };
    }

}


std::tuple<HighPrecision, HighPrecision, HighPrecision> Simulation::runStiffPhase()
{
        /* Find equal time starting from the stiff matter phase. */
        EnergyDensity rhoPhiStiff = phi->energyDensityStiff();
        EnergyDensity rhoStiff = stiff.energyDensity();
    
        auto stiffSolver = EqualTimeSolver(rhoStiff, rhoPhiStiff, p.t0);
        auto [t_eq, rhoStiffEq, rhoPhiEq] = stiffSolver.getEqualTime();

        return std::make_tuple(t_eq, rhoStiffEq, rhoPhiEq);
}

std::tuple<HighPrecision, HighPrecision, HighPrecision> Simulation::runMatterPhase(HighPrecision t0)
{
    EnergyDensity rhoPhiMat = phi->energyDensityMatter(t0);
    EnergyDensity rhoChiMat = chi.energyDensityMatter(t0);
    auto radMatSolver = EqualTimeSolver(rhoPhiMat, rhoChiMat, t0);

    auto [tau_eq, rhoPhiMatEq, rhoChiMatEq] = radMatSolver.getEqualTime();

    return std::make_tuple(tau_eq, rhoPhiMatEq, rhoChiMatEq);
}

std::tuple<HighPrecision, HighPrecision, HighPrecision> Simulation::runRadiationPhase(HighPrecision t0)
{
    EnergyDensity rhoPhiRad = phi->energyDensityRadiation(t0);
    EnergyDensity rhoChiRad = chi.energyDensityRadiation(t0);
    auto radSolver = EqualTimeSolver(rhoPhiRad, rhoChiRad, t0);
    auto [t_eq_rad, rhoPhiRadEq, rhoChiRadEq] = radSolver.getEqualTime();

    return std::make_tuple(t_eq_rad, rhoPhiRadEq, rhoChiRadEq);
}

std::pair<HighPrecision, HighPrecision> Simulation::getReheatingTemperatureAndTime(HighPrecision tau_eq, HighPrecision tau2_eq)
{
    EnergyDensity rhoChiRad = this->chi.energyDensityRadiation(tau_eq);
    auto t_rh = maximize(rhoChiRad, tau2_eq, tau2_eq * HighPrecision("1e5"));
    HighPrecision reheatingTemperature = IntegrationUtils::integrate(this->chi.energyDensityRadiation(tau2_eq), tau2_eq, t_rh);
    HighPrecision T_RH = pow(reheatingTemperature, HighPrecision(1.0/4.0));
    return std::pair(T_RH, t_rh);
}
