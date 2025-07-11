/*
 * ===============================================================================================
 * Project: Reheating via gravitational particle production in the kination epoch
 * Author: J. M. Lankinen
 *
 * Description:
 * This project simulates the decay and evolution of a massive scalar field φ and a massless
 * conformally coupled scalar field χ in the early universe under various cosmological epochs:
 * stiff matter, matter, and radiation domination.
 * It calculates energy densities, equal times, and reheating properties in these epochs.
 * 
 * Purpose:
 * To explore how different parameter choices affect reheating dynamics when inflation ends
 * in kination phase and particles are created gravitationally. Useful for studies in theoretical
 * physics and early universe cosmology.
 *
 * Usage:
 * This file (`main.cpp`) initializes the parameter grid, launches the simulation
 * manager, and coordinates multi-threaded simulation runs. The results are written to a CSV file.
 * The filename can be set in the CSWWriter constructor. You can change the parameter grid in the
 * code below.
 * ===============================================================================================
 */

#include <iostream>
#include <vector>
#include <memory>
#include <chrono>

#include "parameters/parameters.hpp"
#include "simulation/simulation_manager.hpp"
#include "writers/csv_writer.hpp"

using namespace std::chrono;

/*
TODO:
    2. Add support for multiple channels.
*/

int main()
{
    std::vector<ModelParameters> params;
    ModelParameters p;

    std::vector<HighPrecision> lambdaValues{HighPrecision(0.001), HighPrecision(0.01)};
    std::vector<HighPrecision> bValues{HighPrecision(0.01), HighPrecision(0.1),
                                       HighPrecision(1.0), HighPrecision(10),
                                       HighPrecision(100)};
    std::vector<HighPrecision> xiValues{HighPrecision(1.0 / 6.0), HighPrecision(0.0)};
    std::vector<HighPrecision> mValues{};

    // Generate mass points.
    auto generateMassPoints = [&](HighPrecision start, HighPrecision end, int samples)
    {
        HighPrecision step = pow(HighPrecision(10), HighPrecision(1.0) / samples);
        for (HighPrecision m = start; m < end; m*=step)
        {
            mValues.push_back(m);
        }
    };
    // Generate more points in the low mass range where things are interesting.
    generateMassPoints(HighPrecision("1e2"), HighPrecision("1e7"), 100);
    generateMassPoints(HighPrecision("1e7"), HighPrecision("1e25"), 50);


    for (const auto& lambda : lambdaValues)
    {
        for (const auto& xi : xiValues)
        {
            for (const auto& b : bValues)
            {
                for (const auto& m : mValues)
                {
                    p.lambda = lambda;
                    p.b = b;
                    p.xi = xi;
                    p.m = m;
                    params.push_back(p);
                }   
            }
        }
    }

    // Insert filename you want to save results in the parentheses below.
    auto outputWriter = std::make_unique<CSVWriter>("results_dense.csv");

    auto start = steady_clock::now();
    std::cout << "Beginning simulation with " << params.size() << " parameter combinations." << std::endl;

    SimulationManager manager(std::move(params), std::move(outputWriter));
    manager.run(); 
    
    auto end = steady_clock::now();

    auto elapsed = duration_cast<seconds>(end - start);
    auto formatted = hh_mm_ss{elapsed};

    std::cout << "Elapsed time: "
              << formatted.hours().count() << "h "
              << formatted.minutes().count() << "m "
              << formatted.seconds().count() << "s\n";
    
}