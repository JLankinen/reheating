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
     * 4. Add logging
     * 5. Create simulation object. Multithread to each core.
     * 6. Add output to parquet
     * 7. Generate parameter space.
     */

    auto start = std::chrono::high_resolution_clock::now();

    ModelParameters p;
    p.logParameters();

    Simulation sim{p};
    sim.run();

    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::seconds>(end-start).count();
    
    std::cout << "Time taken for one run: " << duration << " seconds." << std::endl;
    
}