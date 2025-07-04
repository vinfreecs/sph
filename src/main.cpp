#include<iostream>
#include"utils/dec.hpp"
#include"utils/particle.hpp"
#include"utils/sph.hpp"
#include"utils/io.hpp"

int main(){
    std::cout<<"Starting SPH sim.."<<std::endl;
    particle p;
    float radius = 0.001;
    float volume = (4.0/3.0) * 3.14 * radius * radius * radius;
    // p.pos_x = 0.0;p.pos_y = 0.0;p.pos_z = 0.0;
    // p.vel_x = 0.0;p.vel_y = 0.0;p.vel_z = 0.0;
    // p.mass = 1.0;
    // p.density = p.mass/volume;
    // p.pressure = 0.0;
    // p.force_x = 0.0;p.force_y = 0.0;p.force_z = 0.0;
    int nParticles =100;
    std::vector<particle> particles(nParticles);
    readInput("../input/pos_100_1.txt",nParticles,particles);
    // for(auto p : particles){
    //     std::cout<<"positions"<<std::endl;
    //     std::cout<<p.pos_x<<" "<<p.pos_y<< " "<<p.pos_z<<std::endl;
    //     std::cout<<"velocities"<<std::endl;
    //     std::cout<<p.vel_x<<" "<<p.vel_y<< " "<<p.vel_z<<std::endl;
    // }
    sph_float dt =0.01;
    sph_float time = 2;
    for(int t=0;t<time/dt;t++){
        calculate_density(particles);
        calculate_force(particles);
        update_velocity(particles,dt);
        update_positions(particles,dt);
        std::string file_name = generateVTKFilename(t+1);
        file_name = "../output/"+file_name;
        writeVTK(file_name, particles);
    }
}