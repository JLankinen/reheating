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
    double reheating_temp;
    double reheating_time;
    double t_eq;
    double rhoStiff_t_eq;
    double rhoPhiStiff_t_eq;
    double rhoChi_t_eq;
    double tau_eq;
    double rhoPhiMatEq;
    double rhoChiMatEq;
    bool toMatter;
    bool bothFound;
};


class Simulation
{
    private:
        const ModelParameters& p;
        std::shared_ptr<PhiParticle> phi;
        ChiParticle chi;
        StiffMatter stiff;
        bool toMatter(const EnergyDensity &rhoChi, double rhoPhi, double timeEquality);
        std::tuple<double, double, double> runMatterPhase(double t0);
        std::tuple<double, double, double> runRadiationPhase(double t0);
        std::pair<double, double> getReheatingTemperatureAndTime(double t_eq);
        std::tuple<bool, double, double, double, bool> runStiffPhase();

    public:
        Simulation(const ModelParameters& p_);
        SimulationResults run();
};


#endif