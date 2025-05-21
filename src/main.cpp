#include <boost/multiprecision/cpp_dec_float.hpp>
#include <iostream>
#include <iomanip>
#include <functional>
#include <chrono>
#include "parameters/parameters.hpp"
#include "model/creation_decay.hpp"
#include "model/energy_densities.hpp"

using HighPrecision = boost::multiprecision::cpp_dec_float_100;

// For simple testing on time.
template <typename F, typename ...Args>
auto TimeCounterDecorator(F&& f, Args&&... args)
{
    auto start = std::chrono::high_resolution_clock::now();
    auto result = std::invoke(std::forward<F>(f), std::forward<Args>(args)...);
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    std::cout << "Elapsed time: " << duration.count()  << " µs" << std::endl;
    return result;
}


int main()
{

    ModelParameters p;
    p.xi = 1.0 / 6.0;
    p.m = HighPrecision("1e20");
    p.t0 = HighPrecision("1e-32");
    p.b = 1;
    p.lambda = 0.01;

    HighPrecision t = HighPrecision("1e-30");
    auto res = RhoChiStiff(p, t);

    std::cout << res;
}