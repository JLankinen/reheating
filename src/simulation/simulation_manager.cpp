#include "simulation/simulation_manager.hpp"

SimulationManager::SimulationManager(std::vector<ModelParameters> params,
                                      std::unique_ptr<ResultsWriter> writer_,
                                      std::size_t workerCount)
    : tasks(std::make_move_iterator(params.begin()),
            std::make_move_iterator(params.end())),
      writer{std::move(writer_)},
      totalSimulationCount{static_cast<int>(params.size())}
    {
        if (workerCount == 0)
        {
            workerCount = 1;
        }
        workers.reserve(workerCount);
    }; 


void SimulationManager::run()
{
    // Launch threads
    for (std::size_t i = 0; i < workers.capacity(); i++)
    {
        workers.emplace_back(&SimulationManager::workerLoop, this);
    }

    for (auto &th : workers)
    {
        th.join();  // Wait until finished
    }
}

void SimulationManager::workerLoop()
{
    while (true)
    {
        ModelParameters p;

        {
            std::unique_lock<std::mutex> lock(queueMtx);
            // Continue if predicate true. Run until simulations are empty or stop condition is met.
            cv.wait(lock, [&] {return stop || !tasks.empty();});

            if (tasks.empty())
            {
                if (stop) {return;}
                else {continue;}
            }

            p = std::move(tasks.front());
            tasks.pop_front();

            // If last task
            if (tasks.empty())
            {
                stop = true;  // Only last thread sets stop.
                cv.notify_all();
            }

        }

        try
        {
            Simulation sim(std::move(p));
            SimulationResults res = sim.run(); // Results of one individual run.
            writer->write(res); // Append the result file.

            int currentSim = ++simulationCounter;
            {
                static std::mutex outputMtx;
                std::lock_guard<std::mutex> lock(outputMtx);
                std::cout << "Simulation: " << currentSim << "/" << totalSimulationCount << std::endl;
            }
        }
        catch (const boost::wrapexcept<std::domain_error>& ex)
        {
            int currentSim = ++simulationCounter;
            std::cerr << "Domain error in simulation: " << ex.what() << "\n";
        }
        catch (const std::exception& ex)
        {
            int currentSim = ++simulationCounter;
            std::cerr << "Standard exception: " << ex.what() << "\n";
        }
        catch (...)
        {
            int currentSim = ++simulationCounter;
            std::cerr << "Unknown error occurred during simulation.\n";
        }
    }
}