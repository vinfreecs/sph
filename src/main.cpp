#include <iostream>
#include <chrono>
#include<fstream>
#include<vector>
#include<iomanip>
#include <fstream>
#include <iostream>
#include <vector>
#include <cmath>
#include "utils/config.hpp"
#include "utils/particle.hpp"
#include "utils/io.hpp"
#include "utils/sph.hpp"


int main()
{
    std::cout << "Starting SPH hash sim.." << std::endl;
    // write_particles_cube_format("../input/pos_150.txt", 150, 0.5);
    int nParticles = 1000;
    std::vector<particle> particles(nParticles);
    readInput("../input/pos_10.txt", nParticles, particles);
    int boundary = 3;
    int nBuckets = ceil(boundary / r_e);
    std::vector<std::vector<int>> buckets(nBuckets * nBuckets * nBuckets);

    std::vector<int> neighbours(particles.size(), -1);
    std::vector<int> particle_neighbours(particles.size());

    
   // update_neighbors(particles, neighbours, particle_neighbours);
    
    // for (int i = 0; i < neighbours.size(); i++)
    // {
    //     std::cout << neighbours[i] << " ";

    // }
    // std::cout<<std::endl;
    // for (int i = 0; i < neighbours.size(); i++)
    // {
    //     std::cout << particle_neighbours[i] << " ";

    // }
    // std::cout<<std::endl;


    sph_float dt = 0.01;
    sph_float time = 10;
    std::string file_name_out = generateVTKFilename(0);
    file_name_out = "../output/" + file_name_out;
    writeVTK(file_name_out, particles);
    sph_float total_time = 0;
    for (int t = 0; t < time / dt; t++)
    {
        auto start = std::chrono::high_resolution_clock::now();
        update_neighbors(particles, neighbours, particle_neighbours);
        
        calculation_density_pressure_hash(particles, neighbours, particle_neighbours);
        // for(auto p:particles){
        //     std::cout<<" Id: "<<p.id<<" Density : "<<p.density<<" pressure : "<<p.pressure<<" Force : "<<p.force_x<<std::endl;
        // }
        reset_forces(particles);
        calculate_force_hash(particles, neighbours, particle_neighbours);
        update_velocity(particles, dt);
        update_positions(particles, dt);
        auto end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed = end - start;
        total_time += elapsed.count();
        std::string file_name = generateVTKFilename(t + 1);
        file_name = "../output/" + file_name;
        writeVTK(file_name, particles);
    }
    std::cout << "Total time taken for simulation : " << total_time << " seconds" << std::endl;
    std::cout << "Simulation completed." << std::endl;
}