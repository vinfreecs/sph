#include <fstream>
#include <iostream>
#include <vector>
#include "particle.hpp"


void write_particles_cube_format(
    const std::string& filename,
    int particles_per_axis,
    double cube_size,
    double mass = 1.0,
    double vx = 0.0, double vy = 0.0, double vz = 0.0
) {
    std::ofstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error: Could not open file " << filename << "\n";
        return;
    }

    double spacing = cube_size / (particles_per_axis - 1);
    int total_particles = particles_per_axis * particles_per_axis * particles_per_axis;

    // Write positions
    file << "positions\n";
    for (int i = 0; i < particles_per_axis; ++i) {
        for (int j = 0; j < particles_per_axis; ++j) {
            for (int k = 0; k < particles_per_axis; ++k) {
                double x = i * spacing;
                double y = j * spacing;
                double z = k * spacing;
                file << x << " " << y << " " << z << "\n";
            }
        }
    }

    // Write masses
    file << "masses\n";
    for (int i = 0; i < total_particles; ++i) {
        file << mass << "\n";
    }

    // Write velocities
    file << "velocities\n";
    for (int i = 0; i < total_particles; ++i) {
        file << vx << " " << vy << " " << vz << "\n";
    }

    file.close();
    std::cout << "Wrote " << total_particles << " particles to " << filename << "\n";
}
