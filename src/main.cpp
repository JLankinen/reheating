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
    std::vector<ModelParameters> params;
    ModelParameters p;

    std::vector<HighPrecision> bValues{HighPrecision(1.0)};
    std::vector<HighPrecision> mValues{HighPrecision("1e2")};

    for (const auto& b : bValues)
    {
        p.b = b;
     for (const auto& m : mValues)
     {
        p.m = m;
        params.push_back(p);
     }   
    }

    auto start = std::chrono::high_resolution_clock::now();
    auto outputWriter = std::make_unique<CSVWriter>("results.csv");

    SimulationManager manager(std::move(params), std::move(outputWriter));
    manager.run();

    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::seconds>(end-start).count();
    
    std::cout << "Time taken for one run: " << duration << " seconds." << std::endl;
    
}