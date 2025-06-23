#include <iostream>
#include <vector>
#include <memory>

#include "parameters/parameters.hpp"
#include "simulation/simulation_manager.hpp"
#include "writers/csv_writer.hpp"


int main()
{
    /**TODO
     * 2. Add logging
     * 3. Generate parameter space.
     */
    std::vector<ModelParameters> params;
    ModelParameters p;

    std::vector<HighPrecision> bValues{HighPrecision(1.0)};
    std::vector<HighPrecision> xiValues{HighPrecision(1.0 / 6.0), HighPrecision(0.0)};
    std::vector<HighPrecision> mValues{HighPrecision("1e2"),
                                       HighPrecision("1e5")};

    for (const auto& xi : xiValues)
    {
        for (const auto& b : bValues)
        {
            for (const auto& m : mValues)
            {
                p.b = b;
                p.xi = xi;
                p.m = m;
                params.push_back(p);
            }   
        }
    }

    auto outputWriter = std::make_unique<CSVWriter>("results.csv");

    SimulationManager manager(std::move(params), std::move(outputWriter));
    manager.run();
    
}