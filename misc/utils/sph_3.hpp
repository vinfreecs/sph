#pragma once
#include "particle.hpp"
#include "spatial_hashing.hpp"
#include "dec.hpp"
#include <algorithm>
#include <cmath>

// Optimized distance calculation - avoid sqrt when possible
inline sph_float distance_squared(const particle& p, const particle& q) {
    const sph_float dx = p.pos_x - q.pos_x;
    const sph_float dy = p.pos_y - q.pos_y;
    const sph_float dz = p.pos_z - q.pos_z;
    return dx*dx + dy*dy + dz*dz;
}

// Precomputed constants for better performance
static const sph_float r_e_squared = r_e * r_e;
static const sph_float r_e_6 = std::pow(r_e, 6);
static const sph_float r_e_9 = std::pow(r_e, 9);
static const sph_float pressure_coeff = 45.0 / (PI * r_e_6);
static const sph_float viscosity_coeff = 45.0 / (PI * r_e_6);
static const sph_float density_coeff = 315.0 / (64.0 * PI * r_e_9);

// Optimized smoothing functions using precomputed constants
inline sph_float smoothing_pressure(sph_float r) {
    if (r < r_e) {
        const sph_float diff = r_e - r;
        return pressure_coeff * diff * diff * diff / r;
    }
    return 0.0;
}

inline sph_float smoothing_viscosity(sph_float r) {
    if (r < r_e) {
        return viscosity_coeff * (r_e - r);
    }
    return 0.0;
}

// Fixed density calculation - was taking sqrt twice
inline sph_float smoothing_density(sph_float r_squared) {
    if (r_squared < r_e_squared) {
        const sph_float diff = r_e_squared - r_squared;
        return density_coeff * diff * diff * diff;
    }
    return 0.0;
}

// More efficient boundary handling
inline void boundary_optimized(particle& p, sph_float boundary = 3.0, sph_float bounce_factor = 0.25) {
    const sph_float half_boundary = boundary * 0.5;
    const sph_float min_pos = -half_boundary + radius;
    const sph_float max_pos = half_boundary - radius;

    // X-axis
    if (p.pos_x < min_pos) {
        p.pos_x = min_pos;
        p.vel_x *= -bounce_factor;
    } else if (p.pos_x > max_pos) {
        p.pos_x = max_pos;
        p.vel_x *= -bounce_factor;
    }

    // Y-axis
    if (p.pos_y < min_pos) {
        p.pos_y = min_pos;
        p.vel_y *= -bounce_factor;
    } else if (p.pos_y > max_pos) {
        p.pos_y = max_pos;
        p.vel_y *= -bounce_factor;
    }

    // Z-axis
    if (p.pos_z < min_pos) {
        p.pos_z = min_pos;
        p.vel_z *= -bounce_factor;
    } else if (p.pos_z > max_pos) {
        p.pos_z = max_pos;
        p.vel_z *= -bounce_factor;
    }
}

// Optimized force calculation
inline void pressure_viscosity_force_optimized(particle& p, particle& q) {
    const sph_float r_squared = distance_squared(p, q);
    
    if (r_squared >= r_e_squared) return; // Early exit for distant particles
    
    const sph_float r = std::sqrt(r_squared);
    const sph_float r_safe = std::max(r, 1e-5f); // Avoid division by zero
    
    const sph_float w_r_pressure = smoothing_pressure(r_safe);
    const sph_float w_r_viscosity = smoothing_viscosity(r_safe);
    
    // Direction vector
    const sph_float dx = q.pos_x - p.pos_x;
    const sph_float dy = q.pos_y - p.pos_y;
    const sph_float dz = q.pos_z - p.pos_z;
    
    // Pressure force calculation (corrected sign)
    const sph_float pressure_magnitude = q.mass * (p.pressure + q.pressure) * w_r_pressure / (2.0 * q.density);
    const sph_float force_pressure_x = pressure_magnitude * dx;
    const sph_float force_pressure_y = pressure_magnitude * dy;
    const sph_float force_pressure_z = pressure_magnitude * dz;
    
    // Viscosity force calculation
    const sph_float viscosity_magnitude = q.mass * w_r_viscosity / q.density;
    const sph_float force_viscosity_x = viscosity_magnitude * (q.vel_x - p.vel_x);
    const sph_float force_viscosity_y = viscosity_magnitude * (q.vel_y - p.vel_y);
    const sph_float force_viscosity_z = viscosity_magnitude * (q.vel_z - p.vel_z);
    
    // Total force
    const sph_float fx = -force_pressure_x + dynamic_viscosity * force_viscosity_x;
    const sph_float fy = -force_pressure_y + dynamic_viscosity * force_viscosity_y;
    const sph_float fz = -force_pressure_z + dynamic_viscosity * force_viscosity_z;
    
    // Apply forces (Newton's third law)
    p.force_x += fx;
    p.force_y += fy;
    p.force_z += fz;
    q.force_x -= fx;
    q.force_y -= fy;
    q.force_z -= fz;
}

// Optimized velocity update
inline void update_velocity_optimized(std::vector<particle>& particles, sph_float dt) {
    for (auto& p : particles) {
        const sph_float inv_density = 1.0 / p.density;
        p.vel_x += p.force_x * dt * inv_density;
        p.vel_y += p.force_y * dt * inv_density;
        p.vel_z += p.force_z * dt * inv_density;
    }
}

// Optimized position update
inline void update_positions_optimized(std::vector<particle>& particles, sph_float dt) {
    for (auto& p : particles) {
        p.pos_x += p.vel_x * dt;
        p.pos_y += p.vel_y * dt;
        p.pos_z += p.vel_z * dt;
        boundary_optimized(p);
    }
}

// Optimized pressure calculation
inline void calculate_pressure_optimized(std::vector<particle>& particles, 
                                       sph_float k = 3.0, 
                                       sph_float rho_o = 998.29, 
                                       sph_float rest_pressure = 10.0) {
    for (auto& p : particles) {
        p.pressure = rest_pressure + k * (p.density - rho_o);
    }
}

// Improved neighbor hash calculation with better cache locality
inline std::vector<int> get_neighbours_hash_optimized(const particle& p, int table_size) {
    const sph_float offset = boundary * 0.5;
    const sph_float inv_r_e = 1.0 / r_e;
    
    const int x = static_cast<int>((p.pos_x + offset) * inv_r_e);
    const int y = static_cast<int>((p.pos_y + offset) * inv_r_e);
    const int z = static_cast<int>((p.pos_z + offset) * inv_r_e);

    std::vector<int> neighbor_hashes;
    neighbor_hashes.reserve(27); // Maximum possible neighbors in 3x3x3 grid

    for (int dx = -1; dx <= 1; dx++) {
        for (int dy = -1; dy <= 1; dy++) {
            for (int dz = -1; dz <= 1; dz++) {
                const int nx = x + dx;
                const int ny = y + dy;
                const int nz = z + dz;

                if (nx >= 0 && ny >= 0 && nz >= 0) {
                    int hash = get_hash(nx, ny, nz);
                    hash = ((hash % table_size) + table_size) % table_size;
                    neighbor_hashes.push_back(hash);
                }
            }
        }
    }

    return neighbor_hashes;
}

// Optimized density calculation using spatial hashing
void calculate_density_hash_optimized(std::vector<particle>& particles,
                                    const std::vector<Neighbour>& neighbours,
                                    const std::vector<int>& start_idx) {
    const int n = particles.size();

    // Reset densities
    for (auto& p : particles) {
        p.density = 0.0;
    }

    // Calculate densities using spatial hashing
    for (int i = 0; i < n; ++i) {
        particle& p = particles[i];
        const std::vector<int> neighbor_cells = get_neighbours_hash_optimized(p, n);

        for (int cell_idx : neighbor_cells) {
            const int start = start_idx[cell_idx];
            if (start == -1) continue;

            for (int k = start; k < neighbours.size() && neighbours[k].cell_id == cell_idx; ++k) {
                const particle& q = particles[neighbours[k].particle_id];
                const sph_float r2 = distance_squared(p, q);
                const sph_float w = smoothing_density(r2);
                p.density += q.mass * w;
            }
        }
    }
}

// Optimized force calculation using spatial hashing
void calculate_force_hash_optimized(std::vector<particle>& particles,
                                  const std::vector<Neighbour>& neighbours,
                                  const std::vector<int>& start_idx) {
    const int n = particles.size();

    // Reset forces
    for (auto& p : particles) {
        p.force_x = p.force_y = p.force_z = 0.0;
    }

    // Calculate inter-particle forces
    for (int i = 0; i < n; ++i) {
        particle& p = particles[i];
        const std::vector<int> neighbor_cells = get_neighbours_hash_optimized(p, n);

        for (int cell_idx : neighbor_cells) {
            const int start = start_idx[cell_idx];
            if (start == -1) continue;

            for (int k = start; k < neighbours.size() && neighbours[k].cell_id == cell_idx; ++k) {
                const int q_idx = neighbours[k].particle_id;
                if (i < q_idx) { // Avoid double computation
                    particle& q = particles[q_idx];
                    pressure_viscosity_force_optimized(p, q);
                }
            }
        }

        // Apply gravity
        p.force_z -= p.mass * gravity;
    }
}

// Improved spatial hashing with better hash distribution
struct Neighbour_Optimized {
    int particle_id = -1;
    int cell_id = -1;

    Neighbour_Optimized() = default;
    Neighbour_Optimized(int pid, int cid) : particle_id(pid), cell_id(cid) {}

    bool operator<(const Neighbour_Optimized& other) const {
        return cell_id < other.cell_id;
    }
};

// Improved hash function with better distribution
inline int get_hash_optimized(int cell_x, int cell_y, int cell_z) {
    // Use larger primes for better distribution
    const int p1 = 15485863;
    const int p2 = 22801763;
    const int p3 = 1291394886;
    
    return (cell_x * p1) ^ (cell_y * p2) ^ (cell_z * p3);
}

// Optimized neighbor update with better memory layout
inline void update_neighbors_optimized(std::vector<particle>& particles,
                                     std::vector<Neighbour_Optimized>& neighbours,
                                     std::vector<int>& start_idx) {
    const int n = particles.size();
    neighbours.resize(n);

    const sph_float offset = boundary * 0.5;
    const sph_float inv_r_e = 1.0 / r_e;

    // Compute cell IDs for all particles
    for (int i = 0; i < n; ++i) {
        const particle& p = particles[i];
        const int x = static_cast<int>((p.pos_x + offset) * inv_r_e);
        const int y = static_cast<int>((p.pos_y + offset) * inv_r_e);
        const int z = static_cast<int>((p.pos_z + offset) * inv_r_e);
        
        const int hash = get_hash_optimized(x, y, z);
        const int cell_id = std::abs(hash % n);
        neighbours[i] = Neighbour_Optimized(i, cell_id);
    }

    // Sort neighbors by cell ID
    std::sort(neighbours.begin(), neighbours.end());

    // Build start index lookup table
    start_idx.clear();
    start_idx.resize(n, -1);

    for (int i = 0; i < n; ++i) {
        const int cell_id = neighbours[i].cell_id;
        if (i == 0 || cell_id != neighbours[i - 1].cell_id) {
            start_idx[cell_id] = i;
        }
    }
}

