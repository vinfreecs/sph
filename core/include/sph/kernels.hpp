#pragma once
// Smoothing kernels (Müller et al. 2003, as used by Harada et al. 2007).

#include "common.hpp"

namespace sph {

// Poly6 kernel W(r), evaluated from r^2:  315 / (64 pi h^9) * (h^2 - r^2)^3
SPH_HD inline sph_float kernel_poly6(sph_float r2, sph_float h) {
    sph_float h2 = h * h;
    if (r2 >= h2) return 0;
    sph_float d = h2 - r2;
    sph_float h3 = h * h * h;
    return 315.0f * d * d * d / (64.0f * kPi * h3 * h3 * h3);
}

// Magnitude of the spiky kernel gradient:  45 / (pi h^6) * (h - r)^3
SPH_HD inline sph_float kernel_spiky_grad(sph_float r, sph_float h) {
    if (r >= h) return 0;
    sph_float d = h - r;
    sph_float h3 = h * h * h;
    return 45.0f * d * d * d / (kPi * h3 * h3);
}

// Viscosity kernel Laplacian:  45 / (pi h^6) * (h - r)
SPH_HD inline sph_float kernel_visc_lap(sph_float r, sph_float h) {
    if (r >= h) return 0;
    sph_float h3 = h * h * h;
    return 45.0f * (h - r) / (kPi * h3 * h3);
}

}  // namespace sph
