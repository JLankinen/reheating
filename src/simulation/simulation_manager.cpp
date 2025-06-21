#include "simulation/simulation_manager.hpp"

SimulationManager::SimulationManager(std::vector<ModelParameters> params,
                                      std::unique_ptr<ResultsWriter> writer_,
                                      std::size_t workerCount)
    : tasks(std::make_move_iterator(params.begin()),
            std::make_move_iterator(params.end())),
      writer(std::move(writer_))
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
            cv.wait(lock, [&] {return stop || !tasks.empty();}); // Run until simulations ar empty.

            if (tasks.empty())
            {
                if (stop) {return;}
                else {continue;}
            }

            p = std::move(tasks.front());
            tasks.pop_front();

            if (tasks.empty())
            {
                stop = true;  // Only last thread sets stop.
                cv.notify_all();
            }

        }

        Simulation sim(p);
        SimulationResults res = sim.run(); // Results of one individual run.
        writer->write(res); // Append the result file.
    }
}