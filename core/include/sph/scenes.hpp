#pragma once
// Particle generators for fluid blocks and ghost-particle boundaries (host only).

#include <cmath>
#include <vector>

#include "particle.hpp"

namespace sph {

inline Particle make_boundary(sph_float x, sph_float y, sph_float z, sph_float spacing_radius) {
    Particle p;
    p.pos_x = x;
    p.pos_y = y;
    p.pos_z = z;
    p.mass = 1.0f;
    p.density = 7870.0f;
    p.fixed = true;
    p.radius = spacing_radius / 2.0f;
    return p;
}

// n fluid particles on a cubic lattice of edge `size`, centred in x/y and sitting above z = size/2.
inline std::vector<Particle> fluid_cube_3d(int n, double size, sph_float radius, sph_float mass = 1.0f) {
    std::vector<Particle> out;
    int per_axis = static_cast<int>(std::ceil(std::cbrt(static_cast<double>(n))));
    double s = size / per_axis;
    for (int i = 0; i < per_axis && (int)out.size() < n; ++i)
        for (int j = 0; j < per_axis && (int)out.size() < n; ++j)
            for (int k = 0; k < per_axis && (int)out.size() < n; ++k) {
                Particle p;
                p.pos_x = sph_float(i * s + s / 2 - size / 2);
                p.pos_y = sph_float(j * s + s / 2 - size / 2);
                p.pos_z = sph_float(k * s + s / 2 + size / 2);
                p.mass = mass;
                p.radius = radius;
                out.push_back(p);
            }
    return out;
}

// Hollow cube of ghost particles with edge `size`, centred at the origin.
inline std::vector<Particle> boundary_box_3d(sph_float size, sph_float spacing_radius, sph_float overlap = 0.2f) {
    std::vector<Particle> out;
    sph_float spacing = spacing_radius * overlap;
    int m = static_cast<int>(size / spacing);
    sph_float half = size / 2.0f;
    for (int i = 0; i <= m; ++i)
        for (int j = 0; j <= m; ++j)
            for (int k = 0; k <= m; ++k) {
                bool border = i == 0 || i == m || j == 0 || j == m || k == 0 || k == m;
                if (border)
                    out.push_back(make_boundary(-half + i * spacing, -half + j * spacing, -half + k * spacing,
                                                spacing_radius));
            }
    return out;
}

// n fluid particles on a square lattice of edge `size`, centred at the origin.
inline std::vector<Particle> fluid_square_2d(int n, double size, sph_float radius, sph_float mass = 1.0f) {
    std::vector<Particle> out;
    int per_axis = static_cast<int>(std::ceil(std::sqrt(static_cast<double>(n))));
    double s = size / per_axis;
    for (int i = 0; i < per_axis && (int)out.size() < n; ++i)
        for (int j = 0; j < per_axis && (int)out.size() < n; ++j) {
            Particle p;
            p.pos_x = sph_float(i * s + s / 2 - size / 2);
            p.pos_y = sph_float(j * s + s / 2 - size / 2);
            p.mass = mass;
            p.radius = radius;
            out.push_back(p);
        }
    return out;
}

// Square outline of ghost particles with edge `size`, centred at the origin.
inline std::vector<Particle> boundary_box_2d(sph_float size, sph_float spacing_radius, sph_float overlap = 0.2f) {
    std::vector<Particle> out;
    sph_float spacing = spacing_radius * overlap;
    int m = static_cast<int>(size / spacing);
    sph_float half = size / 2.0f;
    for (int i = 0; i <= m; ++i)
        for (int j = 0; j <= m; ++j)
            if (i == 0 || i == m || j == 0 || j == m)
                out.push_back(make_boundary(-half + i * spacing, -half + j * spacing, 0, spacing_radius));
    return out;
}

// Ring of ghost particles with an optional gap of `gap_angle` radians starting at `gap_start`.
inline std::vector<Particle> boundary_ring_2d(sph_float ring_radius, sph_float spacing_radius,
                                              sph_float overlap = 0.2f, sph_float gap_angle = 0,
                                              sph_float gap_start = 0) {
    std::vector<Particle> out;
    sph_float spacing = spacing_radius * overlap;
    int count = static_cast<int>(2.0f * kPi * ring_radius / spacing);
    sph_float step = 2.0f * kPi / count;
    for (int i = 0; i < count; ++i) {
        sph_float a = i * step;
        sph_float rel = std::fmod(a - gap_start + 2 * kPi, 2 * kPi);
        if (rel < gap_angle) continue;
        out.push_back(make_boundary(ring_radius * std::cos(a), ring_radius * std::sin(a), 0, spacing_radius));
    }
    return out;
}

inline void append(std::vector<Particle>& dst, const std::vector<Particle>& src) {
    dst.insert(dst.end(), src.begin(), src.end());
}

inline void assign_ids(std::vector<Particle>& ps) {
    for (size_t i = 0; i < ps.size(); ++i) ps[i].id = static_cast<int>(i);
}

}  // namespace sph
