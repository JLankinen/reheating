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

    std::vector<double> lambdaValues{0.1, 0.01, 0.001, 0.0001, 0.00001, 0.000001, 0.0000001};
    std::vector<double> bValues{10.0, 1.0, 0.1};
    std::vector<double> xiValues{0.0, 1.0 / 6.0};  // Minimal and conformal coupling
    std::vector<double> mValues{};


    // Generate mass points.
    auto generateMassPoints = [&](double start, double end, int samples)
    {
        double step = pow(10.0, 1.0 / samples);
        for (double m = start; m < end; m*=step)
        {
            mValues.push_back(m);
        }
    };

    // Generate more points in the low mass range where things are interesting.
    generateMassPoints(1e0, 1e9, 100);  // Low range
    generateMassPoints(1e9, 1e28, 60);  // High range


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
    auto outputWriter = std::make_unique<CSVWriter>("test.csv");

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