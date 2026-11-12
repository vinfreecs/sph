#pragma once
// Per-particle SPH operations shared by the serial CPU loop and the CUDA kernels.
//
// Neighbour search is a hashed cell-linked list with one table slot per particle:
//   head[slot] = last particle inserted into that slot (or -1)
//   next[i]    = previous particle in the same slot      (or -1)
// The CPU inserts sequentially; CUDA inserts with atomicExch (see backends/cuda).
//
// D is the spatial dimension (2 or 3). 2D scenes simply ignore z.

#include "common.hpp"
#include "kernels.hpp"
#include "params.hpp"
#include "particle.hpp"

namespace sph {

SPH_HD inline int cell_coord(sph_float x, sph_float h) { return static_cast<int>(floorf(x / h)); }

SPH_HD inline int cell_slot(int cx, int cy, int cz, int table_size) {
    unsigned hash = (unsigned)cx * 9283711u ^ (unsigned)cy * 689287499u ^ (unsigned)cz * 283923481u;
    return static_cast<int>(hash % (unsigned)table_size);
}

template <int D>
SPH_HD inline int home_slot(const Particle& p, sph_float h, int table_size) {
    int cz = (D == 3) ? cell_coord(p.pos_z, h) : 0;
    return cell_slot(cell_coord(p.pos_x, h), cell_coord(p.pos_y, h), cz, table_size);
}

template <int D>
SPH_HD inline sph_float distance2(const Particle& a, const Particle& b) {
    sph_float dx = a.pos_x - b.pos_x;
    sph_float dy = a.pos_y - b.pos_y;
    sph_float dz = (D == 3) ? a.pos_z - b.pos_z : 0;
    return dx * dx + dy * dy + dz * dz;
}

// Calls visit(j) for every particle j stored in the 3^D cells around p (including p itself).
// Neighbouring cells that hash to the same slot are visited only once.
template <int D, class Visitor>
SPH_HD inline void for_each_neighbour(const Particle& p, const int* head, const int* next, int n,
                                      sph_float h, Visitor&& visit) {
    constexpr int kCells = (D == 3) ? 27 : 9;
    int seen[kCells];
    int n_seen = 0;
    int cx = cell_coord(p.pos_x, h), cy = cell_coord(p.pos_y, h);
    int cz = (D == 3) ? cell_coord(p.pos_z, h) : 0;
    const int zr = (D == 3) ? 1 : 0;
    for (int dx = -1; dx <= 1; ++dx)
        for (int dy = -1; dy <= 1; ++dy)
            for (int dz = -zr; dz <= zr; ++dz) {
                int slot = cell_slot(cx + dx, cy + dy, cz + dz, n);
                bool dup = false;
                for (int s = 0; s < n_seen; ++s) dup |= (seen[s] == slot);
                if (dup) continue;
                seen[n_seen++] = slot;
                for (int j = head[slot]; j != -1; j = next[j]) visit(j);
            }
}

// Density by summation, then pressure from the linear equation of state p = p0 + k (rho - rho0).
template <int D>
SPH_HD inline void compute_density_pressure(int i, Particle* ps, const int* head, const int* next, int n,
                                            const Params& P) {
    Particle& p = ps[i];
    if (p.fixed) {
        p.density = P.boundary_density;
    } else {
        sph_float rho = 0;
        for_each_neighbour<D>(p, head, next, n, P.h, [&](int j) {
            rho += ps[j].mass * kernel_poly6(distance2<D>(p, ps[j]), P.h);
        });
        p.density = rho;
    }
    p.pressure = P.p0 + P.k * (p.density - P.rho0);
}

// Hybrid boundary: if fluid particle p overlaps ghost particle q, push p out along the
// contact normal and remove the normal component of its velocity.
template <int D>
SPH_HD inline void resolve_penetration(Particle& p, const Particle& q) {
    sph_float min_dist = p.radius + q.radius;
    sph_float dist = sqrtf(distance2<D>(p, q));
    if (dist >= min_dist || dist <= 0) return;
    sph_float inv = 1.0f / dist;
    sph_float nx = (p.pos_x - q.pos_x) * inv;
    sph_float ny = (p.pos_y - q.pos_y) * inv;
    sph_float nz = (D == 3) ? (p.pos_z - q.pos_z) * inv : 0;
    sph_float depth = min_dist - dist;
    p.pos_x += nx * depth;
    p.pos_y += ny * depth;
    p.pos_z += nz * depth;
    sph_float vn = p.vel_x * nx + p.vel_y * ny + p.vel_z * nz;
    p.vel_x -= vn * nx;
    p.vel_y -= vn * ny;
    p.vel_z -= vn * nz;
}

// Pressure + viscosity force of q acting on p (accumulated into p only).
template <int D>
SPH_HD inline void add_pair_force(Particle& p, const Particle& q, const Params& P) {
    sph_float r = sqrtf(distance2<D>(p, q));
    if (r > P.h) return;
    if (r < 1e-5f) r = 1e-5f;
    if (q.fixed && P.hybrid_boundary) resolve_penetration<D>(p, q);

    sph_float dx = q.pos_x - p.pos_x;
    sph_float dy = q.pos_y - p.pos_y;
    sph_float dz = (D == 3) ? q.pos_z - p.pos_z : 0;

    sph_float press = q.mass * (p.pressure + q.pressure) * kernel_spiky_grad(r, P.h) / (2.0f * q.density * r);
    sph_float visc = P.viscosity * q.mass * kernel_visc_lap(r, P.h) / q.density;

    p.force_x += -press * dx + visc * (q.vel_x - p.vel_x);
    p.force_y += -press * dy + visc * (q.vel_y - p.vel_y);
    p.force_z += -press * dz + visc * (q.vel_z - p.vel_z);
}

template <int D>
SPH_HD inline void compute_forces(int i, Particle* ps, const int* head, const int* next, int n,
                                  const Params& P) {
    Particle& p = ps[i];
    if (p.fixed) return;
    p.force_x = p.force_y = p.force_z = 0;
    for_each_neighbour<D>(p, head, next, n, P.h, [&](int j) {
        if (j != i) add_pair_force<D>(p, ps[j], P);
    });
    sph_float fg = p.density * P.gravity;
    if (P.gravity_axis == 0) p.force_x -= fg;
    else if (P.gravity_axis == 1) p.force_y -= fg;
    else p.force_z -= fg;
}

SPH_HD inline void clamp_axis(sph_float& x, sph_float& v, sph_float r, sph_float half, sph_float bounce) {
    if (x - r < -half) { x = -half + r; v *= -bounce; }
    else if (x + r > half) { x = half - r; v *= -bounce; }
}

// Explicit (symplectic) Euler step.
template <int D>
SPH_HD inline void integrate(Particle& p, sph_float dt, const Params& P) {
    if (p.fixed) return;
    if (p.density > 1e-5f) {
        p.vel_x += p.force_x * dt / p.density;
        p.vel_y += p.force_y * dt / p.density;
        p.vel_z += p.force_z * dt / p.density;
    }
    p.pos_x += p.vel_x * dt;
    p.pos_y += p.vel_y * dt;
    p.pos_z += p.vel_z * dt;
    if (P.wall_half > 0) {
        sph_float r = P.particle_radius;
        clamp_axis(p.pos_x, p.vel_x, r, P.wall_half, P.wall_bounce);
        clamp_axis(p.pos_y, p.vel_y, r, P.wall_half, P.wall_bounce);
        if (D == 3) clamp_axis(p.pos_z, p.vel_z, r, P.wall_half, P.wall_bounce);
    }
}

}  // namespace sph
