#ifndef RESULTS_WRITER_H_
#define RESULTS_WRITER_H_

#include "simulation/simulation.hpp"

/**
 * @brief Abstract interface for writing simulation results.
 *
 * Defines a generic interface for outputting SimulationResults.
 */
class ResultsWriter
{
    public:
        virtual ~ResultsWriter() = default;

        /**
         * @brief Writes a single simulation result to the output.
         * 
         * @param res The result data of a simulation to be written.
         */
        virtual void write(const SimulationResults& res) = 0;
};

#endif