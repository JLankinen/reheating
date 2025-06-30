# Description
This project simulates the decay and evolution of a massive scalar field φ and a massless
conformally coupled scalar field χ in the early universe under various cosmological epochs:
stiff matter, matter, and radiation domination.
It calculates energy densities, equal times, and reheating properties in these epochs.
This is a C++ implementation of the Python code in https://arxiv.org/abs/1910.07520.

# Installation
Use CMake to install. Needs boost library and GoogleTest

```bash
mkdir build
cd build
cmake ../CMakeLists.txt
make
```