#include<fstream>
#include<vector>
#include<iomanip>
#include<iostream>
#include"particle.hpp"



void writeParticlesToVTK(const std::string &filename, std::vector<particle> particles, int nParticles)
{
    std::ofstream vtkfile(filename);
    if (!vtkfile)
    {
        std::cerr << "Error writing to file: " << filename << "\n";
        return;
    }

    vtkfile << "# vtk DataFile Version 4.0\n";
    vtkfile << "Simple particle visualization\n";
    vtkfile << "ASCII\n";
    vtkfile << "DATASET UNSTRUCTURED_GRID\n";

    vtkfile << "POINTS " << nParticles << " double\n";
    for (int i = 0; i < nParticles; ++i)
    {
        vtkfile << std::fixed << std::setprecision(6)
                << particles[i].pos_x << " " << particles[i].pos_y << " " << particles[i].pos_z << "\n";
    }

    vtkfile << "CELLS 0 0\n";
    vtkfile << "CELL_TYPES 0\n";

    vtkfile << "POINT_DATA " << nParticles << "\n";

    vtkfile << "SCALARS m double 1\n";
    vtkfile << "LOOKUP_TABLE default\n";
    for (int i = 0; i < nParticles; ++i)
    {
        vtkfile << particles[i].mass << "\n";
    }

    vtkfile << "VECTORS v double\n";
    for (int i = 0; i < nParticles; ++i)
    {
        vtkfile << particles[i].vel_x << " " << particles[i].vel_y << " " << particles[i].vel_z << "\n";
    }

    vtkfile.close();
}

std::string generateVTKFilename(int frame)
{
    return "frame_" + std::to_string(frame) + ".vtk";
}

void writeVTK(const std::string &filename, std::vector<particle> &particles)
{
    std::ofstream out(filename);
    if (!out.is_open())
    {
        std::cerr << "Failed to open file: " << filename << std::endl;
        return;
    }

    const size_t numParticles = particles.size();

    out << "# vtk DataFile Version 4.0\n";
    out << "hesp visualization file\n";
    out << "ASCII\n";
    out << "DATASET UNSTRUCTURED_GRID\n";
    out << "POINTS " << numParticles << " double\n";

    out << std::fixed << std::setprecision(6);
    for (auto particle : particles)
    {
        double x = particle.pos_x;
        double y = particle.pos_y;
        double z = particle.pos_z;
        out << x << " " << y << " " << z << "\n";
    }

    out << "CELLS 0 0\n";
    out << "CELL_TYPES 0\n";
    out << "POINT_DATA " << numParticles << "\n";

    // Write masses
    out << "SCALARS m double\n";
    out << "LOOKUP_TABLE default\n";
    for (auto m : particles)
    {
        out << m.mass << "\n";
    }

    // Write velocities
    out << "VECTORS v double\n";
    for (auto particle : particles)
    {
        double vx = particle.vel_x;
        double vy = particle.vel_y;
        double vz = particle.vel_z;
        out << vx << " " << vy << " " << vz << "\n";
    }

    out.close();
}

void readInput(const std::string &filename, int nParticles, std::vector<particle> &particles)
{
    std::ifstream file(filename);
    std::string line;

    if (!file)
    {
        std::cerr << "Error opening file: " << filename << "\n";
        return;
    }

    // Read until "positions"
    while (std::getline(file, line))
    {
        if (line.find("positions") != std::string::npos)
            break;
    }

    for (int i = 0; i < nParticles; ++i)
    {
        file >> particles[i].pos_x >> particles[i].pos_y >> particles[i].pos_z;
    }

    // Read until "masses"
    while (std::getline(file, line))
    {
        if (line.find("masses") != std::string::npos)
            break;
    }

    for (int i = 0; i < nParticles; ++i)
    {
        file >> particles[i].mass;
    }

    // Read until "velocities"
    while (std::getline(file, line))
    {
        if (line.find("velocities") != std::string::npos)
            break;
    }

    for (int i = 0; i < nParticles; ++i)
    {
        file >> particles[i].vel_x >> particles[i].vel_y >> particles[i].vel_z;

    }

    for(int i=0; i<nParticles;i++){
        particles[i].id = i;
    }
    

    file.close();
}