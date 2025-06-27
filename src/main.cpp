#include <iostream>
#include <vector>
#include <memory>
#include <chrono>

#include "parameters/parameters.hpp"
#include "simulation/simulation_manager.hpp"
#include "writers/csv_writer.hpp"

using namespace std::chrono;

int main()
{
    /**TODO
     * 1. Add failed parameters file.
     * 2. Add logging
     */
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
    generateMassPoints(HighPrecision("1e1"), HighPrecision("1e8"), 100);
    generateMassPoints(HighPrecision("1e8"), HighPrecision("1e25"), 40);


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

    auto outputWriter = std::make_unique<CSVWriter>("results_dense.csv");

    std::cout << "Beginning simulation with " << params.size() << " parameter combinations.";
    auto start = steady_clock::now();

    SimulationManager manager(std::move(params), std::move(outputWriter));
    //manager.run(); 
    
    auto end = steady_clock::now();

    auto elapsed = duration_cast<seconds>(end - start);
    auto formatted = hh_mm_ss{elapsed};

    std::cout << "Elapsed time: "
              << formatted.hours().count() << "h "
              << formatted.minutes().count() << "m "
              << formatted.seconds().count() << "s\n";
    
}