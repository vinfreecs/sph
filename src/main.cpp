#include <iostream>
#include <chrono>
#include <fstream>
#include <vector>
#include <iomanip>
#include <fstream>
#include <iostream>
#include <vector>
#include <cmath>
#include "utils/config_serial.hpp"
#include "utils/particle.hpp"
#include "utils/io.hpp"
#include "utils/sph.hpp"

int main()
{
    std::cout << "Starting SPH hash sim.." << std::endl;
    write_particles_cube_format("../input/pos_10000.txt", 10000, 0.5);
    int nParticles = 10000;
    std::vector<particle> particles(nParticles);
    readInput("../input/pos_10000.txt", nParticles, particles);

    std::vector<particle> boundary_vec = create_boundary_particles_3d(1.5, r_e_b);
    particles.insert(particles.end(), boundary_vec.begin(), boundary_vec.end());

    std::vector<int> neighbours(particles.size(), -1);
    std::vector<int> particle_neighbours(particles.size(),-1);


    sph_float dt = 0.001;
    sph_float time = 5;
    std::string file_name_out = generateVTKFilename(0);
    file_name_out = "../output/" + file_name_out;
    writeVTK(file_name_out, particles);
    sph_float total_time = 0;
    update_neighbors(particles, neighbours, particle_neighbours);

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
        if (t % 200 == 0)
        {

            std::string file_name = generateVTKFilename(t + 1);
            file_name = "../output/" + file_name;

            auto start_write_vtk = std::chrono::high_resolution_clock::now();
            writeVTK(file_name, particles);
            auto end_write_vtk = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> elapsed_write_vtk = end_write_vtk - start_write_vtk;
            sph_float total_time_write_vtk = elapsed_write_vtk.count();
            std::cout << "Time Taken to write vtk :" << total_time_write_vtk << std::endl;
            std::cout << "Time Taken for "<<t <<" steps : "<< total_time << std::endl;

        }
    }
    std::cout << "Total time taken for simulation : " << total_time << " seconds" << std::endl;
    std::cout << "Simulation completed." << std::endl;
}