#pragma once
#include <vector>
#include"dec.hpp"

struct particle
{
    int id;
    float x, y;
    float radius;
    float mass;
    float pressure;
    float density;
    float velocity_x, velocity_y;
    int bucket_id = -1;
    bool is_boundary = false;
    // particle type fluid and boundary
};

std::vector<particle> create_boundary_particles(
    sph_float box_width,
    sph_float box_height,
    sph_float spacing,
    sph_float radius,
    sph_float mass)
{
    std::vector<particle> boundaries;

    spacing = 2.0f * radius;  // ensure tight packing BEFORE row/col calculation

    int rows = static_cast<int>(box_height / spacing) + 1;
    int cols = static_cast<int>(box_width / spacing) + 1;
    float offset_x = -0.5f * box_width;
    float offset_y = -0.5f * box_height;

    // Bottom and Top boundaries
    for (int j = 0; j < cols; ++j)
    {
        for (int i : {0, rows - 1})
        {
            particle p;
            p.x = offset_x + j * spacing;
            p.y = offset_y + i * spacing;
            p.radius = radius;
            p.mass = mass;
            p.density = 0.0f;
            p.pressure = 0.0f;
            p.velocity_x = 0.0f;
            p.velocity_y = 0.0f;
            p.bucket_id = -1;
            p.is_boundary = true;
            p.id = i * cols + j;
            boundaries.push_back(p);
        }
    }

    // Left and Right boundaries (excluding corners)
    for (int i = 1; i < rows - 1; ++i)
    {
        for (int j : {0, cols - 1})
        {
            particle p;
            p.x = offset_x + j * spacing;
            p.y = offset_y + i * spacing;
            p.radius = radius;
            p.mass = mass;
            p.density = 0.0f;
            p.pressure = 0.0f;
            p.velocity_x = 0.0f;
            p.velocity_y = 0.0f;
            p.bucket_id = -1;
            p.is_boundary = true;
            p.id = i * cols + j;
            boundaries.push_back(p);
        }
    }

    return boundaries;
}

std::vector<particle> create_particles_grid(int rows, int cols, float spacing, float radius, float mass)
{
    std::vector<particle> particles;
    float offset_x = -0.5f * (cols - 1) * spacing;
    float offset_y = -0.5f * (rows - 1) * spacing;
    spacing = 2.0f*radius;

    for (int i = 0; i < rows; ++i)
    {
        for (int j = 0; j < cols; ++j)
        {
            particle p;
            p.x = offset_x + j * spacing;
            p.y = offset_y + i * spacing;
            p.radius = radius;
            p.mass = mass;
            p.density = 0.0f;
            p.pressure = 0.0f;
            p.velocity_x = 0.0f;
            p.velocity_y = 0.0f;
            p.id = i * rows + j;
            p.bucket_id = -1;
            particles.push_back(p);
        }
    }
    return particles;
}

std::vector<particle> create_all_particles()
{
    std::vector<particle> fluid_particles = create_particles_grid(
        20, 20,
        0.05f,
        0.01f,
        1.0f);

    std::vector<particle> boundary_particles = create_boundary_particles(
        2.0f, 2.0f, // width, height of the container box
        0.018f,      // spacing
        0.002f,      // radius
        1.0f        // mass
    );
    std::vector<particle> particles;
    particles.insert(particles.end(), fluid_particles.begin(), fluid_particles.end());
    particles.insert(particles.end(), boundary_particles.begin(), boundary_particles.end());

    return particles;
}