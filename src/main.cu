#include <iostream>
#include <chrono>
#include <fstream>
#include <vector>
#include <iomanip>
#include <cmath>
#include <string>
#include <cuda_runtime.h>
#include "utils/config.hpp"
#include "utils/particle.hpp"
#include "utils/io.hpp"
#include "utils/sph_cuda.hpp"
#include "utils/cuda_error.hpp"


int main()
{
    std::cout << "Starting SPH hash sim.." << std::endl;
    int nParticles = 50000;
    std::string filename = "../input/pos_"+std::to_string(nParticles)+".txt";
    write_particles_cube_format(filename, nParticles, 0.5); //to write fluid particles in a cube format
    std::vector<particle> particles(nParticles);
    readInput(filename, nParticles, particles);

    //to get boundary particles and attach them to particles
    std::vector<particle> boundary_vec = create_boundary_particles_3d(1.5, r_e_b);
    particles.insert(particles.end(), boundary_vec.begin(), boundary_vec.end());

    //initializing neighors/cell head and next
    std::vector<int> neighbours(particles.size(), -1);
    std::vector<int> particle_neighbours(particles.size(),-1);

    nParticles = particles.size();

    //intitializing array for cuda /device
    particle *d_particles;
    int *d_neighbours;
    int *d_particle_neighbours;

    size_t size_particles = nParticles * sizeof(particle);
    size_t size_int = nParticles * sizeof(int);

    //allocating memory and copying the dat from host to device
    checkCudaError(cudaMalloc(&d_particles, size_particles));
    checkCudaError(cudaMalloc(&d_neighbours, size_int));
    checkCudaError(cudaMalloc(&d_particle_neighbours, size_int));

    checkCudaError(cudaMemcpy(d_particles, particles.data(), size_particles, cudaMemcpyHostToDevice));
    checkCudaError(cudaMemcpy(d_neighbours, neighbours.data(), size_int, cudaMemcpyHostToDevice));
    checkCudaError(cudaMemcpy(d_particle_neighbours, particle_neighbours.data(), size_int, cudaMemcpyHostToDevice));


    
    int threadsPerBlock = 256;
    int blocksPerGrid = (nParticles + threadsPerBlock - 1) / threadsPerBlock;

    //the time step and total time
    sph_float dt = time_step;
    sph_float time = total_time;
    
    //writing the first orentation to get an idea of the orentation
    std::string file_name_out = generateVTKFilename(0);
    file_name_out = "../output/" + file_name_out;
    writeVTK(file_name_out, particles);
    sph_float total_time = 0;
    for (int t = 0; t < time / dt; t++)
    {
        auto start = std::chrono::high_resolution_clock::now();

        //cudaMemset to initialise the neighbors/cell head and next with -1 or else can write a kernal to intialize to -1
        checkCudaError(cudaMemset(d_neighbours, -1, size_int));
        checkCudaError(cudaMemset(d_particle_neighbours, -1, size_int));

        //updating the cells 
        update_neighbors_kernel<<<blocksPerGrid, threadsPerBlock>>>(d_particles, d_neighbours, d_particle_neighbours, nParticles);
        checkCudaError(cudaDeviceSynchronize());

        //kernal to calculate density and then the pressure
        compute_density_pressure_hash_kernel<<<blocksPerGrid, threadsPerBlock>>>(d_particles, d_neighbours, d_particle_neighbours, nParticles);
        checkCudaError(cudaDeviceSynchronize());
        // removed reset as alredy setting force to zero in force kernel
        //calculating internal(pressure and viscous) and external forces(gravity)
        calculate_force_hash_kernel<<<blocksPerGrid, threadsPerBlock>>>(d_particles, d_neighbours, d_particle_neighbours, nParticles);
        checkCudaError(cudaDeviceSynchronize());

        //updating the velocity and position using euler method
        update_velocity_kernel<<<blocksPerGrid, threadsPerBlock>>>(d_particles, nParticles, dt);
        checkCudaError(cudaDeviceSynchronize());
        update_positions_kernel<<<blocksPerGrid, threadsPerBlock>>>(d_particles, nParticles, dt);
        checkCudaError(cudaDeviceSynchronize());
        auto end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed = end - start;
        total_time += elapsed.count();


        //to write vtk output for every x number of iterations
        if (t % output_freq == 0)
        {
            cudaMemcpy(particles.data(), d_particles, size_particles, cudaMemcpyDeviceToHost);
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

    //imp to free the allocated mem
    checkCudaError(cudaFree(d_particles));
    checkCudaError(cudaFree(d_neighbours));
    checkCudaError(cudaFree(d_particle_neighbours));
    std::cout << "Total time taken for simulation : " << total_time << " seconds" << std::endl;
    std::cout << "Simulation completed." << std::endl;
}