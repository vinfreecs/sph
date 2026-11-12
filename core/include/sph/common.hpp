#pragma once
// Shared definitions for every backend (CPU, CUDA, OpenGL, Metal).
// Everything in core/ must compile both as plain C++17 and under nvcc.

#include <cmath>

#if defined(__CUDACC__)
#define SPH_HD __host__ __device__
#else
#define SPH_HD
#endif

namespace sph {

using sph_float = float;  // switch to double here to change precision everywhere

constexpr sph_float kPi = 3.14159265358979323846f;

}  // namespace sph
