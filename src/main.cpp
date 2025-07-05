#include<iostream>
#include "utils/dec.hpp"
#include"utils/dec.hpp"
#include"utils/particle.hpp"
#include"utils/sph.hpp"
#include"utils/io.hpp"
#include"utils/create_particle.hpp"

int main(){
    std::cout<<"Starting SPH sim.."<<std::endl;
    particle p;
    sph_float volume = (4.0/3.0) * 3.14 * radius * radius * radius;
    // p.pos_x = 0.0;p.pos_y = 0.0;p.pos_z = 0.0;
    // p.vel_x = 0.0;p.vel_y = 0.0;p.vel_z = 0.0;
    // p.mass = 1.0;
    // p.density = p.mass/volume;
    // p.pressure = 0.0;
    // p.force_x = 0.0;p.force_y = 0.0;p.force_z = 0.0;

    write_particles_cube_format("../input/pos_8.txt", 2, 1);
    int nParticles =1000;
    std::vector<particle> particles(nParticles);
    readInput("../input/pos_10.txt",nParticles,particles);
    // for(auto p : particles){
    //     std::cout<<"positions"<<std::endl;
    //     std::cout<<p.pos_x<<" "<<p.pos_y<< " "<<p.pos_z<<std::endl;
    //     std::cout<<"velocities"<<std::endl;
    //     std::cout<<p.vel_x<<" "<<p.vel_y<< " "<<p.vel_z<<std::endl;
    // }
    int boundary = 3;
    int nBuckets = boundary/r_e;
    //std::cout<<"Number of buckets : "<<nBuckets*nBuckets*nBuckets<<std::endl;
    std::vector<std::vector<int>> buckets(nBuckets*nBuckets*nBuckets); //neighbours 3 is the boundary and 0.006 is the smoothing length effective radius etc
    update_buckets(particles,buckets,boundary,nBuckets);
    sph_float dt =0.01;
    sph_float time = 2;
    std::string file_name_out = generateVTKFilename(0);
    file_name_out = "../output/"+file_name_out;
    writeVTK(file_name_out, particles);
    for(int t=0;t<time/dt;t++){
        update_buckets(particles,buckets,boundary,nBuckets);
        calculate_density(particles);
        calculate_pressure(particles);
        //calculate_force(particles);
        calculate_force_buckets(particles,buckets, nBuckets);
        update_velocity(particles,dt);
        update_positions(particles,dt);
        std::cout<<"Time step : "<<t+1<<std::endl;
        // for(auto p:particles){
        //     std::cout<<"Particle "<<p.id<<" position : "<<p.pos_x<<" "<<p.pos_y<<" "<<p.pos_z<<std::endl;
        //     std::cout<<"Particle "<<p.id<<" velocity : "<<p.vel_x<<" "<<p.vel_y<<" "<<p.vel_z<<std::endl;
        //     std::cout<<"Particle "<<p.id<<" density : "<<p.density<<std::endl;
        //     std::cout<<"Particle "<<p.id<<" pressure : "<<p.pressure<<std::endl;
        // }
        std::string file_name = generateVTKFilename(t+1);
        file_name = "../output/"+file_name;
        writeVTK(file_name, particles);
    }
}