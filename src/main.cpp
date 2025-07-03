#include<iostream>
#include"utils/sph.hpp"

int main(){
    std::cout<<"Starting SPH sim.."<<std::endl;
    particle p;
    float radius = 0.001;
    float volume = (4.0/3.0) * 3.14 * radius * radius * radius;
    p.pos_x = 0.0;p.pos_y = 0.0;p.pos_z = 0.0;
    p.vel_x = 0.0;p.vel_y = 0.0;p.vel_z = 0.0;
    p.mass = 1.0;
    p.density = p.mass/volume;
    p.pressure = 0.0;
    p.force_x = 0.0;p.force_y = 0.0;p.force_z = 0.0;
    std::vector<particle> particles;
    particles.push_back(p);
}