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
#include "utils/io_aks.hpp"
#include "utils/sph_cuda.hpp"
#include "utils/cuda_error.hpp"


int main()
{
    std::cout << "Starting SPH hash sim.." << std::endl;
    // write_particles_cube_format("../input/pos_150.txt", 150, 0.5);
    int nParticles = 94392;
    std::vector<particle> particles(nParticles);
    // readInput("../input/pos_100000.txt", nParticles, particles);
    readVTK("../input/particle_scene.vtk", particles);


    std::vector<int> neighbours(particles.size(), -1);
    std::vector<int> particle_neighbours(particles.size(),-1);

    particle *d_particles;

    int *d_neighbours;

    int *d_particle_neighbours;

    size_t size_particles = nParticles * sizeof(particle);
    size_t size_int = nParticles * sizeof(int);

    checkCudaError(cudaMalloc(&d_particles, size_particles));
    checkCudaError(cudaMalloc(&d_neighbours, size_int));
    checkCudaError(cudaMalloc(&d_particle_neighbours, size_int));

    checkCudaError(cudaMemcpy(d_particles, particles.data(), size_particles, cudaMemcpyHostToDevice));
    checkCudaError(cudaMemcpy(d_neighbours, neighbours.data(), size_int, cudaMemcpyHostToDevice));
    checkCudaError(cudaMemcpy(d_particle_neighbours, particle_neighbours.data(), size_int, cudaMemcpyHostToDevice));



    int threadsPerBlock = 256;
    int blocksPerGrid = (nParticles + threadsPerBlock - 1) / threadsPerBlock;

    sph_float dt = 0.001f;
    sph_float time = 2.0f;
    std::string file_name_out = generateVTKFilename(0);
    file_name_out = "../output/out_boundary_92k/" + file_name_out;
    writeVTK(file_name_out, particles);
    sph_float total_time = 0;
    for (int t = 0; t < time / dt; t++)
    {
        auto start = std::chrono::high_resolution_clock::now();

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

        if (t % 50 == 0)
        {
            cudaMemcpy(particles.data(), d_particles, size_particles, cudaMemcpyDeviceToHost);
            std::string file_name = generateVTKFilename(t + 1);
            file_name = "../output/out_boundary_92k/" + file_name;
            writeVTK(file_name, particles);
        }
    }

    checkCudaError(cudaFree(d_particles));
    checkCudaError(cudaFree(d_neighbours));
    checkCudaError(cudaFree(d_particle_neighbours));
    std::cout << "Total time taken for simulation : " << total_time << " seconds" << std::endl;
    std::cout << "Simulation completed." << std::endl;
}