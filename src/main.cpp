#include <iostream>
#include <iomanip>
#include <chrono>
#include <vector>
#include <memory>

#include "parameters/parameters.hpp"
#include "simulation/simulation_manager.hpp"
#include "writers/csv_writer.hpp"

int main()
{
    /**TODO
     * 1. Fix the restriction
     * (1.2. Add refinement e.g., newton to solver builder)
     * 2. Add logging
     * 3. Generate parameter space.
     */

    auto start = std::chrono::high_resolution_clock::now();

    ModelParameters p;
    p.logParameters();

    std::vector<ModelParameters> params;

    params.push_back(p);

    auto outputWriter = std::make_unique<CSVWriter>();

    SimulationManager manager(std::move(params), std::move(outputWriter));
    manager.run();

    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::seconds>(end-start).count();
    
    std::cout << "Time taken for one run: " << duration << " seconds." << std::endl;
    
}