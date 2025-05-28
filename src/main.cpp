#include <iostream>
#include <iomanip>
#include <chrono>
#include "parameters/parameters.hpp"
#include "solvers/equal_time_solver.hpp"
#include "model/particles/phi_particle.hpp"
#include "model/particles/chi_particle.hpp"
#include "model/particles/stiff_matter.hpp"

int main()
{

    ModelParameters p;
    p.xi = 1.0 / 6.0;
    p.m = HighPrecision("1e40");
    p.t0 = HighPrecision("1e-32");
    p.b = 1;
    p.lambda = 0.01;

    auto phiParticle = PhiParticle(p);
    auto stiff = StiffMatter(p);

    EnergyDensity rhoPhiStiff = phiParticle.energyDensityStiff();
    EnergyDensity rhoStiff = stiff.energyDensity();
    auto chiParticle = ChiParticle(p, rhoStiff);
    EnergyDensity rhoChiStiff = chiParticle.energyDensityStiff();

    auto stiffSolver = EqualTimeSolver(rhoStiff, rhoPhiStiff, 1e-32).withRestriction(rhoChiStiff);
    auto [t_eq, rhoStiffEq, rhoPhiEq] = stiffSolver.getEqualTime();

    std::cout << "t_eq: " << t_eq << " RhoPhi at t_eq: " << rhoPhiEq << " RhoStiff at t_eq: " << rhoStiffEq << std::endl;

}