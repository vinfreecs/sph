#pragma once
#include "dec.hpp"
#include "particle.hpp"
#include <vector>
#include <cmath>

const int boundary = 3;

inline int get_hash(int cell_x, int cell_y, int cell_z)
{
    return (cell_x * 9283711) ^ (cell_y * 689287499) ^ (cell_z * 283923481);
}

inline void update_neighbors(std::vector<particle> &particles,
                             std::vector<int> &neighbours,
                             std::vector<int> &particle_neighbours)
{
    int n = particles.size();
    std::fill(neighbours.begin(), neighbours.end(), -1);
    std::fill(particle_neighbours.begin(), particle_neighbours.end(), -1);

    for (int i = 0; i < n; ++i)
    {
        particle &p = particles[i];
        int x = (int)((p.pos_x + (sph_float)boundary / 2.0f) / r_e);
        int y = (int)((p.pos_y + (sph_float)boundary / 2.0f) / r_e);
        int z = (int)((p.pos_z + (sph_float)boundary / 2.0f) / r_e);
        int hash = get_hash(x, y, z);
        int cell_id = ((hash % n) + n) % n; // use a fixed-size table (size = n here)
        int temp = neighbours[cell_id];
        neighbours[cell_id] = i;
        particle_neighbours[i] = temp;
    }
}

void calculation_density_pressure_hash(std::vector<particle> &particles,
                                  const std::vector<int> &neighbours,
                                  const std::vector<int> &particle_neighbours)
{
    int n = particles.size();
    for (int i = 0; i < n; i++)
    {
        particle &p = particles[i];
        p.density = 0;
        int x = (int)((p.pos_x + (sph_float)boundary / 2.0f) / r_e);
        int y = (int)((p.pos_y + (sph_float)boundary / 2.0f) / r_e);
        int z = (int)((p.pos_z + (sph_float)boundary / 2.0f) / r_e);
        for (int dx = -1; dx <= 1; ++dx)
        {
            for (int dy = -1; dy <= 1; ++dy)
            {
                for (int dz = -1; dz <= 1; ++dz)
                {
                    int hx = x + dx;
                    int hy = y + dy;
                    int hz = z + dz;
                    int hash = get_hash(hx, hy, hz);
                    int cell_id = ((hash % n) + n) % n; // same fix as above

                    int j = neighbours[cell_id];
                    while (j != -1)
                    {
                        if (j != i)
                        {
                            particle &q = particles[j];
                            sph_float r = distance_squared(p, q);
                            if (r < r_e)
                            {
                                sph_float w_r = smoothing_density(r);
                                p.density += q.mass * (w_r);
                            }
                        }
                        j = particle_neighbours[j];
                    }
                }
            }
        }
        p.pressure = rest_pressure + k * (p.density - rho_o);
    }
}

void calculate_force_hash(std::vector<particle> &particles,
                            std::vector<int> &neighbours,
                            std::vector<int> &particle_neighbours)
{
    int n = particles.size();
    for (int i = 0; i < n; i++)
    {
        particle &p = particles[i];
        int x = (int)((p.pos_x + (sph_float)boundary / 2.0f) / r_e);
        int y = (int)((p.pos_y + (sph_float)boundary / 2.0f) / r_e);
        int z = (int)((p.pos_z + (sph_float)boundary / 2.0f) / r_e);
        for (int dx = -1; dx <= 1; ++dx)
        {
            for (int dy = -1; dy <= 1; ++dy)
            {
                for (int dz = -1; dz <= 1; ++dz)
                {
                    int hx = x + dx;
                    int hy = y + dy;
                    int hz = z + dz;
                    int hash = get_hash(hx, hy, hz);
                    int cell_id = ((hash % n) + n) % n; // same fix as above

                    int j = neighbours[cell_id];
                    while (j != -1)
                    {
                        if (j != i)
                        {
                            particle &q = particles[j];
                            sph_float r = distance_squared(p, q);
                            if (r < r_e)
                            {
                                pressure_viscosity_force(p,q);
                            }
                        }
                        j = particle_neighbours[j];
                    }
                }
            }
        }
        p.force_z -= p.mass * gravity; // gravity only in the -ve z direction 

    } 
}

void reset_forces(std::vector<particle>& particles){
    for(auto & p : particles){
        p.force_x=0;
        p.force_y=0;
        p.force_z=0;
    }
}