#pragma once
#include <vector>
#include <iostream>
#include "particle.hpp"
#include "dec.hpp"


void update_buckets(std::vector<particle> &particles, std::vector<std::vector<int>> &buckets, int boundary, int nBuckets)
{
    for (int i = 0; i < particles.size(); i++)
    {
        sph_float x = particles[i].x + (sph_float)boundary / 2;
        sph_float y = particles[i].y + (sph_float)boundary / 2;
        if (x < 0 || x > boundary || y < 0 || y > boundary)
        {
            // std::cerr<<"Particle "<<i<<" is out of bounds. Please check the input file."<<std::endl;
            // std::cerr<<"x: "<<x<<" y: "<<y<<" z: "<<" bucket_idx: "<<std::endl;
            // std::cerr<<"nBuckets: "<<nBuckets<<" boundary: "<<boundary<<std::endl;
            continue;
        }
        int bucket_idx = (int)(x / r_e) + int(y / r_e) * nBuckets;

        if (bucket_idx > nBuckets * nBuckets - 1 || bucket_idx < 0)
        {
            // std::cerr<<"Particle "<<i<<" is out of bounds. Please check the input file."<<std::endl;
            // std::cerr<<"x: "<<x<<" y: "<<y<<" z: "<<" bucket_idx: "<<bucket_idx<<std::endl;
            // std::cerr<<"nBuckets: "<<nBuckets<<" boundary: "<<boundary<<std::endl;
            continue;
        }
        if (particles[i].bucket_id == bucket_idx)
        {
            continue; // already in the same bucket
        }
        else
        {
            if (particles[i].bucket_id != -1)
            {
                // Removing from old bucket
                int old_bucket_idx = particles[i].bucket_id;
                for (int j = 0; j < buckets[old_bucket_idx].size(); j++)
                {
                    if (buckets[old_bucket_idx][j] == particles[i].id)
                    {
                        buckets[old_bucket_idx].erase(buckets[old_bucket_idx].begin() + j);
                        break;
                    }
                }
            }
            particles[i].bucket_id = bucket_idx;
            buckets[bucket_idx].push_back(particles[i].id);
        }

        // std::cout<<"Particle "<<i<<" in bucket "<<bucket_idx<<std::endl;
    }
}


int calculate_cell_hash(particle &p)
{
    int x = std::floor(p.x / r_e);
    int y = std::floor(p.x / r_e);

    int x_hash = x * 73856093;
    int y_hash = y * 19349663;
    return x_hash + y_hash;
    // 83,492,791 for z
}

void update_neighbours_list(std::vector<int> neighbour_list, std::vector<particle> particles)
{
    for (int i = 0; i < particles.size(); i++)
    {
        particle &p = particles[i];
        int cell_ind = calculate_cell_hash(p) % (neighbour_list.size());
        // neighbour_list[i] =
    }
}

std::vector<int> get_neighbours(int idx, int n)
{
    int x = idx % n;
    int y = (idx / n);
    std::vector<int> neighbors;
    for (int dx = -1; dx <= 1; dx++)
    {
        for (int dy = -1; dy <= 1; dy++)
        {
            int nx = x + dx;
            int ny = y + dy;
            if (nx >= 0 && ny >= 0 && nx < n && ny < n)
            {
                int neighbor_idx = nx + ny * n;
                neighbors.push_back(neighbor_idx);
            }
        }
    }
    return neighbors;
}
