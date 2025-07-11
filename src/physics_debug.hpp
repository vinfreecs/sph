#pragma once
#include <vector>
#include <iostream>
#include "utils/particle.hpp"
#include "utils/dec.hpp"
#include "utils/neighbour.hpp"

void apply_boundary_conditions(std::vector<particle> &particles, float window_width, float window_height, float bounce = 1.0f)
{
    for (auto &p : particles)
    {
        // Left or right wall
        if (p.x - p.radius < -window_width / 2.0f)
        {
            p.x = -window_width / 2.0f + p.radius;
            p.velocity_x *= bounce;
        }
        else if (p.x + p.radius > window_width / 2.0f)
        {
            p.x = window_width / 2.0f - p.radius;
            p.velocity_x *= bounce;
        }

        // Bottom or top wall
        if (p.y - p.radius < -window_height / 2.0f)
        {
            p.y = -window_height / 2.0f + p.radius;
            p.velocity_y *= bounce;
        }
        else if (p.y + p.radius > window_height / 2.0f)
        {
            p.y = window_height / 2.0f - p.radius;
            p.velocity_y *= bounce;
        }
    }
}


float smoothing_density(sph_float r)
{
    r = std::sqrt(r);
    if (r < r_e)
    {
        sph_float val = 315 * std::pow((r_e * r_e - r * r), 3) / (64 * PI * std::pow(r_e, 9));
        return val;
    }
    else
        return 0.0;
}

float distance_squared(particle p, particle q)
{
    sph_float dx = p.x - q.x;
    sph_float dy = p.y - q.y;
    sph_float r_dist = dx * dx + dy * dy;
    return r_dist;
}

void calculate_density(std::vector<particle> &particles)
{
    for (auto &p : particles)
    {
        if (p.is_boundary)
            continue;
        p.density = 0.0;
        for (auto &q : particles)
        {
            if (!q.is_boundary)
            {
                float r_dist_squared = distance_squared(p, q);
                float w_r = smoothing_density(r_dist_squared);
                p.density += q.mass * (w_r);
            }
            else
            {
                sph_float r_dist = distance_squared(p, q);
                if (ideal_d * ideal_d < r_dist)
                    continue;
                sph_float w_r = smoothing_density(r_dist);
                p.density += q.mass * w_r;
                // std::cout<<"inside the density boundary calcu - "<<w_r<<" density of the moving particle "<<p.density<<" "<<std::endl;
            }
        }
    }
}

void calculate_pressure(std::vector<particle> &particles, sph_float k = 3, sph_float rho_o = 998.8, sph_float rest_pressure = 0.0)
{
    for (auto &p : particles)
    {
        p.pressure = rest_pressure + k * (p.density - rho_o);
        // std::cout<<"Particle "<<p.id<<" density :"<<p.density<<" pressure :"<<rest_pressure + k * (p.density - rho_o)<<std::endl;
    }
}

sph_float smoothing_pressure(sph_float r)
{
    if (r < 1e-5)
        r = 1e-5; // to avoid division by zero
    if (r < r_e)
    {
        sph_float val = 45 * std::pow(r_e - r, 3) / (PI * std::pow(r_e, 6) * r);
        return val;
    }
    else
        return 0.0;
}

sph_float smoothing_viscosity(sph_float r)
{
    if (r_e > r)
    {
        sph_float val = 45 * (r_e - r) / (PI * std::pow(r_e, 6));
        return val;
    }
    else
        return 0.0;
}

void pressure_viscosity_force(particle &p, particle &q, float dt)
{
    sph_float r_squared = distance_squared(p, q);
    sph_float r = std::sqrt(r_squared);
    if (r > r_e)
        return; // no interaction if distance is greater than smoothing length
    // if(r < 1e-5) return; // to avoid division by zero
    //  Calculate the force due to pressure and viscosity

    sph_float w_r_pressure = smoothing_pressure(r);
    sph_float w_r_viscosity = smoothing_viscosity(r);
    sph_float force_pressure_x = q.mass * (p.pressure + q.pressure) * w_r_pressure * (q.x - p.x) / (2 * q.density); // TODO
    sph_float force_pressure_y = q.mass * (p.pressure + q.pressure) * w_r_pressure * (q.y - p.y) / (2 * q.density);
    sph_float force_viscosity_x = q.mass * (q.velocity_x - p.velocity_x) * w_r_viscosity / (q.density);
    sph_float force_viscosity_y = q.mass * (q.velocity_y - p.velocity_y) * w_r_viscosity / (q.density);
    sph_float fx = -force_pressure_x;
    sph_float fy = -force_pressure_y;

    p.velocity_x += ((-force_pressure_x / p.density) + (dynamic_viscosity * force_viscosity_x / p.density)) * dt;
    p.velocity_y += ((-force_pressure_y / p.density) + (dynamic_viscosity * force_viscosity_y / p.density)) * dt;

    q.velocity_x -= ((-force_pressure_x / p.density) + (dynamic_viscosity * force_viscosity_x / p.density)) * dt;
    q.velocity_y -= ((-force_pressure_y / p.density) + (dynamic_viscosity * force_viscosity_y / p.density)) * dt;
}

void calculate_force(std::vector<particle> &particles, float dt)
{
    // Force = gravity_force + pressure_force + viscosity_force
    for (int i = 0; i < particles.size(); i++)
    {
        particle &p = particles[i];
        if (p.is_boundary)
            continue;
        for (int j = 0; j < particles.size(); j++)
        {
            particle &q = particles[j];
            if (!q.is_boundary && i != j)
            {
                pressure_viscosity_force(p, q, dt);
            }

            if (q.is_boundary)
            {
                sph_float r2 = distance_squared(p, q);
                sph_float riw = std::sqrt(r2);
                if (riw < 1e-5f)
                    riw = 1e-5f; // avoid division by zero
                sph_float dx = p.x - q.x;
                sph_float dy = p.y - q.y;
                sph_float nri_x = dx / riw;
                sph_float nri_y = dy / riw;

                sph_float overlap = ideal_d - riw;
                if (riw > ideal_d) continue;


                // the incoming force or velocity is not considered?

                sph_float force_wall_x = p.mass * overlap * (nri_x) / (dt * dt);
                sph_float force_wall_y = p.mass * overlap * (nri_y) / (dt * dt);

                std::cout << "velocity before x:" << p.velocity_x << "  y:" << p.velocity_y << " " << std::endl;

                std::cout << "boundary repulsion force in x:" << (force_wall_x / p.density) * dt << "  y:" << (force_wall_y / p.density) * dt << " " << std::endl;
                //  Apply repulsion
                std::cout << "Positions of q.x and q.y :" << q.x << "  " << q.y << std::endl;
                p.velocity_x += 20000*(force_wall_x / p.density)*dt;
                p.velocity_y += 20000*(force_wall_y / p.density)*dt;

                std::cout << "velocity after x:" << p.velocity_x << "  y:" << p.velocity_y << " " << std::endl;
            }
        }
        p.velocity_y += -0.8 * dt;
    }
}


void calculate_density_buckets(std::vector<particle> &particles, std::vector<std::vector<int>> &buckets, int nBuckets)
{
    for (auto &p : particles)
    {
        if (p.is_boundary == true)
            continue;
        p.density = 0.0;
        int p_bucket_id = p.bucket_id;
        std::vector<int> neighbour_buckets = get_neighbours(p_bucket_id, nBuckets);
        for (int j = 0; j < neighbour_buckets.size(); j++)
        {
            int neighbour_bucket_id = neighbour_buckets[j];
            if (neighbour_bucket_id < 0 || neighbour_bucket_id >= nBuckets * nBuckets)
                continue;
            for (int k = 0; k < buckets[neighbour_bucket_id].size(); k++)
            {
                particle &q = particles[buckets[neighbour_bucket_id][k]];
                if (!q.is_boundary)
                {
                    float r_dist_squared = distance_squared(p, q);
                    float w_r = smoothing_density(r_dist_squared);
                    p.density += q.mass * (w_r);
                }
                else
                {
                    sph_float r_dist = distance_squared(p, q);
                    if (ideal_d * ideal_d < r_dist)
                        continue;
                    sph_float w_r = smoothing_density(r_dist);
                    p.density += q.mass * w_r;
                    // std::cout<<"inside the density boundary calcu - "<<w_r<<" density of the moving particle "<<p.density<<" "<<std::endl;
                }
            }
        }
    }
}

void calculate_force_buckets(std::vector<particle> &particles, std::vector<std::vector<int>> &buckets, int nBuckets, sph_float dt)
{
    for (int i = 0; i < particles.size(); i++)
    {
        particle &p = particles[i];
        if (p.is_boundary)
            continue;

        int p_bucket_id = p.bucket_id;
        std::vector<int> neighbour_buckets = get_neighbours(p_bucket_id, nBuckets);

        for (int neighbor_bucket_id : neighbour_buckets)
        {
            if (neighbor_bucket_id < 0 || neighbor_bucket_id >= nBuckets * nBuckets)
                continue;

            for (int k : buckets[neighbor_bucket_id])
            {
                particle &q = particles[k];

                if (!q.is_boundary && p.id < q.id)
                {
                    pressure_viscosity_force(p, q, dt);
                }
                // if (q.is_boundary)
                // {
                //     std::cout << "entered" << std::endl;
                //     sph_float r2 = distance_squared(p, q);

                //     sph_float riw = std::sqrt(r2);
                //     // if (riw < 1e-5f) riw = 1e-5f; // avoid division by zero

                //     sph_float dx = p.x - q.x;
                //     sph_float dy = p.y - q.y;
                //     sph_float nri_x = dx / riw;
                //     sph_float nri_y = dy / riw;

                //     sph_float overlap = ideal_d - riw;
                //     if (ideal_d < riw)
                //         continue;

                //     sph_float force_wall_x = p.mass * overlap * (nri_x) / (dt * dt);
                //     sph_float force_wall_y = p.mass * overlap * (nri_y) / (dt * dt);

                //     // std::cout<<"velocity before x:"<<p.velocity_x<<"  y:"<<p.velocity_y<<" "<<std::endl;

                //     // std::cout<<"boundary repulsion force in x:"<<(force_wall_x / p.density)*dt<<"  y:"<<(force_wall_y/ p.density)*dt<<" "<<std::endl;
                //     //  Apply repulsion
                //     if (q.x != 0 || q.x != 2)
                //         p.velocity_x += (force_wall_x)*dt;
                //     if (q.y != 0 || q.y != 2)
                //         p.velocity_y += (force_wall_y)*dt;

                //     std::cout << "velocity after x:" << p.velocity_x << "  y:" << p.velocity_y << " " << std::endl;
                // }
            }
        }
    }
}

void update_positions(std::vector<particle> &particles, float dt)
{
    for (auto &p : particles)
    {
        p.x += p.velocity_x * dt;
        p.y += p.velocity_y * dt;
    }
}
