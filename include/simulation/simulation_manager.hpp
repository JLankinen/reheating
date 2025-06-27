#ifndef SIM_MANAGER_H_
#define SIM_MANAGER_H_

#include <deque>
#include <thread>
#include <mutex>
#include <condition_variable>
#include <atomic>
#include <memory>
#include <vector>

#include "simulation/simulation.hpp"
#include "writers/results_writer.hpp"

/**
 * @brief Manages concurrent execution of simulations and result collection.
 *
 * Coordinates the execution of multiple simulations in parallel using a worker thread pool.
 * Distributes simulation tasks from a shared queue and writes the results to a file.
 */
class SimulationManager
{
    private:
        // Task queue
        std::deque<ModelParameters> tasks;
        std::mutex queueMtx;
        std::condition_variable cv;
        // Workers
        std::vector<std::thread> workers;
        std::atomic<bool> stop{false};
        // Writer
        std::unique_ptr<ResultsWriter> writer;
        // Counter
        std::atomic<int> simulationCounter{0};
        int totalSimulationCount = 0;

        void workerLoop();

    public:
        SimulationManager(std::vector<ModelParameters> params,
                          std::unique_ptr<ResultsWriter> writer,
                          std::size_t workerCount = std::thread::hardware_concurrency());
        void run();
};


#endif