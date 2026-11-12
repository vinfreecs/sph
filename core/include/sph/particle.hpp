#pragma once
#include "common.hpp"

namespace sph {

// One struct for fluid and boundary particles. 2D scenes leave the z fields at 0.
// Boundary ("ghost") particles have fixed == true and never move.
struct Particle {
    int id = 0;
    sph_float radius = 0;
    sph_float pos_x = 0, pos_y = 0, pos_z = 0;
    sph_float vel_x = 0, vel_y = 0, vel_z = 0;
    sph_float mass = 1;
    sph_float density = 0;
    sph_float pressure = 0;
    sph_float force_x = 0, force_y = 0, force_z = 0;
    bool fixed = false;
};

}  // namespace sph
