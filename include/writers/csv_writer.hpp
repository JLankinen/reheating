#ifndef CSV_WRITER_H_
#define CSV_WRITER_H_

#include <mutex>
#include <sstream>
#include <fstream>
#include <string>
#include <filesystem>
#include "writers/results_writer.hpp"

using namespace std::filesystem;

/**
 * @brief Write results to a CSV file.
 *
 * Outputs SimulationResults to a CSV file named "results.csv" inside a "results" directory
 * created one level up from where the script is run.
 */
class CSVWriter : public ResultsWriter
{
    private:
        std::ofstream fs;
        std::mutex mtx;
        path outputDir = current_path().parent_path() / "results";
        path outputFile = outputDir / "results.csv";

        static std::string toCSVRow(const SimulationResults& res);

    public:
        explicit CSVWriter();
        void write(const SimulationResults& res) override;
};


#endif