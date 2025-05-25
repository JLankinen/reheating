#include <iostream>
#include <iomanip>
#include <functional>
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
    auto start2 = std::chrono::high_resolution_clock::now();

    auto phiParticle = PhiParticle(p);
    auto stiff = StiffMatter(p);

    EnergyDensity rhoPhiStiff = phiParticle.energyDensityStiff();
    EnergyDensity rhoStiff = stiff.energyDensity();

    auto chiParticle = ChiParticle(p, rhoStiff);
    EnergyDensity rhoChiStiff = chiParticle.energyDensityStiff();

    auto solver = EqualTimeSolver(rhoStiff, rhoPhiStiff, p, 1e-32, 1e-10).withRestriction(rhoChiStiff);
    HighPrecision t_eq_s = solver.getEqualTime();
    auto end2 = std::chrono::high_resolution_clock::now();
    auto duration2 = std::chrono::duration_cast<std::chrono::microseconds>(end2-start2);
    std::cout << "Time elapsed for solver: " << duration2.count() << std::endl;
    std::cout << "Time of equality for solver: " << t_eq_s;

}