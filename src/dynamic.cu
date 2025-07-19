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
    int maxParticles = 30000;
    std::vector<particle> particles; // host buffer preallocated
    std::vector<particle> boundary_vec = create_boundary_particles_3d(1.5, r_e_b);
    int boundary_particles = boundary_vec.size();
    particles.insert(particles.end(), boundary_vec.begin(), boundary_vec.end());


    std::vector<int> neighbours(particles.size()+maxParticles, -1);
    std::vector<int> particle_neighbours(particles.size()+maxParticles,-1);
    std::cout<<"Num of particles before :"<<nParticles<<std::endl;

    int nParticles = particles.size();
    std::cout<<"Num of particles after :"<<nParticles<<std::endl;


size_t size_particles = (nParticles+maxParticles) * sizeof(particle);
size_t size_int = (nParticles+maxParticles) * sizeof(int);

checkCudaError(cudaMalloc(&d_particles, size_particles));
checkCudaError(cudaMalloc(&d_neighbours, size_int));
checkCudaError(cudaMalloc(&d_particle_neighbours, size_int));


    particle *d_particles;

    int *d_neighbours;

    int *d_particle_neighbours;


    checkCudaError(cudaMalloc(&d_particles, size_particles));
    checkCudaError(cudaMalloc(&d_neighbours, size_int));
    checkCudaError(cudaMalloc(&d_particle_neighbours, size_int));

    checkCudaError(cudaMemcpy(d_particles, particles.data(), size_particles, cudaMemcpyHostToDevice));
    checkCudaError(cudaMemcpy(d_neighbours, neighbours.data(), size_int, cudaMemcpyHostToDevice));
    checkCudaError(cudaMemcpy(d_particle_neighbours, particle_neighbours.data(), size_int, cudaMemcpyHostToDevice));



    int threadsPerBlock = 256;
    int blocksPerGrid = (nParticles + threadsPerBlock - 1) / threadsPerBlock;

    sph_float dt = 0.0001f;
    sph_float time = 5.0f;
    std::string file_name_out = generateVTKFilename(0);
    file_name_out = "../output/" + file_name_out;
    writeVTK(file_name_out, particles);
    sph_float total_time = 0;
    for (int t = 0; t < time / dt; t++)
    {
        auto start = std::chrono::high_resolution_clock::now();

        if (nParticles-boundary_particles < maxParticles) {
    particle newP;
    newP.pos_x = -0.5f;
    newP.pos_y = -0.48f;
    newP.pos_z = 0.0f;
    newP.vel_x = 4.0f;
    newP.vel_y = 0.0f;
    newP.vel_z = 0.0f;
    newP.mass = 1.0f;
    newP.density = 0.0f;
    newP.pressure = 0.0f;
    newP.force_x = 0.0f;
    newP.force_y = 0.0f;
    newP.force_z = 0.0f;
    newP.fixed = false;
    newP.radius = 0.03f;

    // Add to host array
    particles[nParticles] = newP;

    // Copy only this one to device
    checkCudaError(cudaMemcpy(&d_particles[nParticles], &particles[nParticles], sizeof(particle), cudaMemcpyHostToDevice));

    nParticles++; // increment active count
}

        checkCudaError(cudaMemset(d_neighbours, -1, size_int));
        checkCudaError(cudaMemset(d_particle_neighbours, -1, size_int));

        update_neighbors_kernel<<<blocksPerGrid, threadsPerBlock>>>(d_particles, d_neighbours, d_particle_neighbours, nParticles);
        checkCudaError(cudaDeviceSynchronize());
        compute_density_pressure_hash_kernel<<<blocksPerGrid, threadsPerBlock>>>(d_particles, d_neighbours, d_particle_neighbours, nParticles);
        checkCudaError(cudaDeviceSynchronize());
        // removed reset as alredy setting force to zero in force kernel

        calculate_force_hash_kernel<<<blocksPerGrid, threadsPerBlock>>>(d_particles, d_neighbours, d_particle_neighbours, nParticles);
        checkCudaError(cudaDeviceSynchronize());
        update_velocity_kernel<<<blocksPerGrid, threadsPerBlock>>>(d_particles, nParticles, dt);
        checkCudaError(cudaDeviceSynchronize());
        update_positions_kernel<<<blocksPerGrid, threadsPerBlock>>>(d_particles, nParticles, dt);
        checkCudaError(cudaDeviceSynchronize());
        auto end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed = end - start;
        total_time += elapsed.count();

        if (t % 500 == 0)
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

    checkCudaError(cudaFree(d_particles));
    checkCudaError(cudaFree(d_neighbours));
    checkCudaError(cudaFree(d_particle_neighbours));
    std::cout << "Total time taken for simulation : " << total_time << " seconds" << std::endl;
    std::cout << "Simulation completed." << std::endl;
}