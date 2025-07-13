#include<iostream>
#include "utils/dec.hpp"
#include"utils/dec.hpp"
#include"utils/particle.hpp"
#include"utils/sph.hpp"
#include"utils/io.hpp"
#include<chrono>

int main(){
    std::cout<<"Starting SPH sim.."<<std::endl;
    particle p;
    sph_float volume = (4.0/3.0) * 3.14 * radius * radius * radius;

    write_particles_cube_format("../input/pos_5.txt", 5, 0.5);
    int nParticles =1000;
    std::vector<particle> particles(nParticles);
    readInput("../input/pos_10.txt",nParticles,particles);
    int boundary = 3;
    int nBuckets = ceil(boundary/r_e);
    std::vector<std::vector<int>> buckets(nBuckets*nBuckets*nBuckets); //neighbours 3 is the boundary and 0.006 is the smoothing length effective radius etc
    update_buckets(particles,buckets,boundary,nBuckets);
    std::vector<int> neighbour_list(particles.size());
    std::vector<int> start_neighbour(particles.size(),-1);
    sph_float dt =0.01;
    sph_float time = 10;
    std::string file_name_out = generateVTKFilename(0);
    file_name_out = "../output/"+file_name_out;
    writeVTK(file_name_out, particles);
    sph_float total_time =0;
    for(int t=0;t<time/dt;t++){
        auto start = std::chrono::high_resolution_clock::now();
        update_buckets(particles,buckets,boundary,nBuckets);

        //calculate_density(particles);
        calculate_density_buckets(particles, buckets, nBuckets);
        for(auto p:particles){
            std::cout<<" Id: "<<p.id<<" Density : "<<p.density<<" pressure : "<<p.pressure<<" Force : "<<p.force_x<<std::endl;
        }
        calculate_pressure(particles);
        //calculate_force(particles);
        calculate_force_buckets(particles,buckets, nBuckets);
        update_velocity(particles,dt);
        update_positions(particles,dt);
        auto end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed = end - start;
        total_time += elapsed.count();
        std::string file_name = generateVTKFilename(t+1);
        file_name = "../output/"+file_name;
        writeVTK(file_name, particles);
    }
    std::cout<<"Total time taken for simulation : "<<total_time<<" seconds"<<std::endl;
    std::cout<<"Simulation completed."<<std::endl;
}