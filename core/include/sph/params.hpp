#pragma once
#include "common.hpp"

namespace sph {

// Physical and numerical parameters. Passed by value to CUDA kernels, so keep it trivially copyable.
struct Params {
    sph_float h = 0.04f;                  // smoothing radius r_e
    sph_float viscosity = 10.0f;          // mu in the viscosity force
    sph_float gravity = 9.8f;
    int gravity_axis = 2;                 // 0 = x, 1 = y, 2 = z (gravity points in -axis)
    sph_float k = 3.0f;                   // pressure stiffness: p = p0 + k (rho - rho0)
    sph_float rho0 = 998.29f;             // rest density
    sph_float p0 = 0.0f;                  // rest pressure
    sph_float boundary_density = 999.0f;  // density assigned to fixed (ghost) particles
    sph_float particle_radius = 0.03f;    // radius given to fluid particles
    sph_float boundary_spacing = 0.02f;   // r_e_b, spacing scale of generated boundary particles
    bool hybrid_boundary = true;          // push fluid out of penetrated ghost particles + reflect velocity
    sph_float wall_half = 0.0f;           // > 0 enables a fallback clamp to the box [-wall_half, wall_half]
    sph_float wall_bounce = 0.25f;        // restitution of that fallback clamp
};

// 3D closed box (HESP submission): 10k fluid particles, ghost + hybrid boundary.
inline Params preset_box_3d() { return Params{}; }

// 3D scenes loaded from VTK (tools/make_scene.py): free fall / flow over sphere runs.
// These runs used gravity along -y, higher viscosity and a non-zero rest pressure.
inline Params preset_vtk_3d() {
    Params p;
    p.viscosity = 40.0f;
    p.gravity_axis = 1;
    p.p0 = 10.0f;
    p.particle_radius = 0.01f;
    return p;
}

// 2D "flow in a cup" (OpenGL demo, slide scene 3).
inline Params preset_cup_2d() {
    Params p;
    p.viscosity = 0.0001f;
    p.gravity_axis = 1;
    p.rho0 = 1000.0f;
    p.particle_radius = 0.008f;
    p.boundary_spacing = 0.008f;
    p.hybrid_boundary = false;
    p.wall_half = 1.0f;
    return p;
}

}  // namespace sph
