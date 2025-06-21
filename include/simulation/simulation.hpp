#ifndef SIMULATION_H_
#define SIMULATION_H_

#include <memory>

#include "model/energy/creation_decay.hpp"
#include "model/particles/chi_particle.hpp"
#include "model/particles/phi_particle.hpp"
#include "model/particles/stiff_matter.hpp"
#include "parameters/parameters.hpp"


struct SimulationResults
{
    // Original model parameters
    ModelParameters params;
    // Calculated quantities
    HighPrecision reheating_temp;
    HighPrecision reheating_time;
    HighPrecision t_eq;
    HighPrecision rhoStiff_t_eq;
    HighPrecision rhoPhiStiff_t_eq;
    HighPrecision tau_eq;
    HighPrecision rhoPhiMatEq;
    HighPrecision rhoChiMatEq;
    HighPrecision tau2_eq;
    HighPrecision rhoPhiRadEq;
    HighPrecision rhoChiRadEq;
};


class Simulation
{
    private:
        const ModelParameters& p;
        std::shared_ptr<PhiParticle> phi;
        ChiParticle chi;
        StiffMatter stiff;
        std::tuple<HighPrecision, HighPrecision, HighPrecision> runStiffPhase();
        std::tuple<HighPrecision, HighPrecision, HighPrecision> runMatterPhase(HighPrecision t0);
        std::tuple<HighPrecision, HighPrecision, HighPrecision> runRadiationPhase(HighPrecision t0);
        std::pair<HighPrecision, HighPrecision> getReheatingTemperatureAndTime(HighPrecision rhoChiRad,
                                                                               HighPrecision t_eq);
    public:
        Simulation(const ModelParameters& p_);
        SimulationResults run();
};


#endif