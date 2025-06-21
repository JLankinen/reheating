#include <ios>

#include "writers/csv_writer.hpp"

#include <iostream>
#include <filesystem>

CSVWriter::CSVWriter()
{
    if(!std::filesystem::exists(outputDir))
    {
        std::filesystem::create_directory(outputDir);
    }

    if (!std::filesystem::exists(outputFile))
    {
        std::ofstream fs(outputFile);
    }

    fs.open(outputFile, std::ios::app);
    if(!fs.is_open())
    {
        throw std::runtime_error("Cannot open csv file.");
    }
    
    // Write the csv header
    fs << "t0,m,lambda,b,xi,G_N,"
            "reheating_temp,reheating_time,"
            "t_eq,rhoStiff_t_eq,rhoPhiStiff_t_eq,"
            "tau_eq,rhoPhiMatEq,rhoChiMatEq,"
            "tau2_eq,rhoPhiRadEq,rhoChiRadEq\n";
};

void CSVWriter::write(const SimulationResults& res)
{
    std::lock_guard<std::mutex> lock(mtx);
    fs << toCSVRow(res) << "\n";
}

std::string CSVWriter::toCSVRow(const SimulationResults& r)
{
    std::ostringstream ss;
    ss << r.params.t0      << ',' << r.params.m     << ',' << r.params.lambda << ','
    << r.params.b       << ',' << r.params.xi    << ',' << r.params.G_N    << ','
    << r.reheating_temp << ',' << r.reheating_time << ','
    << r.t_eq           << ',' << r.rhoStiff_t_eq    << ',' << r.rhoPhiStiff_t_eq << ','
    << r.tau_eq         << ',' << r.rhoPhiMatEq      << ',' << r.rhoChiMatEq      << ','
    << r.tau2_eq        << ',' << r.rhoPhiRadEq      << ',' << r.rhoChiRadEq;
    return ss.str();
}