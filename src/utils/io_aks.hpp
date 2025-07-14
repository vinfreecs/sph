#include <fstream>
#include <sstream>
#include <iostream>
#include <string>
#include <vector>
#include <iomanip>
//#include "utils/particle.hpp"  // Make sure this includes your 'particle' struct

// Generate filename like frame_0.vtk
std::string generateVTKFilename(int frame)
{
    return "frame_" + std::to_string(frame) + ".vtk";
}

// Write particle data to VTK
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
    for (const auto& particle : particles)
    {
        out << particle.pos_x << " "
            << particle.pos_y << " "
            << particle.pos_z << "\n";
    }

    out << "CELLS 0 0\n";
    out << "CELL_TYPES 0\n";
    out << "POINT_DATA " << numParticles << "\n";

    out << "SCALARS m double\nLOOKUP_TABLE default\n";
    for (const auto& particle : particles)
    {
        out << particle.mass << "\n";
    }

    out << "VECTORS v double\n";
    for (const auto& particle : particles)
    {
        out << particle.vel_x << " "
            << particle.vel_y << " "
            << particle.vel_z << "\n";
    }

    out.close();
}

// Read particle data from a VTK file written by your Python script
void readVTK(const std::string &filename, std::vector<particle> &particles)
{
    std::ifstream file(filename);
    std::string line;

    if (!file)
    {
        std::cerr << "Error opening file: " << filename << "\n";
        return;
    }

    int nParticles = 0;

    // === Skip to "POINTS" line ===
    while (std::getline(file, line))
    {
        if (line.find("POINTS") != std::string::npos)
        {
            std::istringstream iss(line);
            std::string keyword, dtype;
            iss >> keyword >> nParticles >> dtype;
            break;
        }
    }

    if (nParticles <= 0)
    {
        std::cerr << "Error: POINTS section not found or invalid.\n";
        return;
    }

    particles.resize(nParticles);

    // === Read positions ===
    for (int i = 0; i < nParticles; ++i)
    {
        if (!(file >> particles[i].pos_x >> particles[i].pos_y >> particles[i].pos_z))
        {
            std::cerr << "Error reading position for particle " << i << "\n";
            return;
        }
        particles[i].id = i;
    }

    // === Skip to "SCALARS mass" ===
    while (std::getline(file, line))
    {
        if (line.find("SCALARS mass") != std::string::npos)
            break;
    }
    std::getline(file, line); // skip LOOKUP_TABLE line

    for (int i = 0; i < nParticles; ++i)
    {
        if (!(file >> particles[i].mass))
        {
            std::cerr << "Error reading mass for particle " << i << "\n";
            return;
        }
    }

    // === Skip to "VECTORS velocity" ===
    while (std::getline(file, line))
    {
        if (line.find("VECTORS velocity") != std::string::npos)
            break;
    }

    for (int i = 0; i < nParticles; ++i)
    {
        if (!(file >> particles[i].vel_x >> particles[i].vel_y >> particles[i].vel_z))
        {
            std::cerr << "Error reading velocity for particle " << i << "\n";
            return;
        }
    }

    // === Skip to "SCALARS radius" ===
    while (std::getline(file, line))
    {
        if (line.find("SCALARS radius") != std::string::npos)
            break;
    }
    std::getline(file, line); // skip LOOKUP_TABLE line

    for (int i = 0; i < nParticles; ++i)
    {
        if (!(file >> particles[i].radius))
        {
            std::cerr << "Error reading radius for particle " << i << "\n";
            return;
        }
    }

    // === Skip to "SCALARS fixed" ===
    while (std::getline(file, line))
    {
        if (line.find("SCALARS fixed") != std::string::npos)
            break;
    }
    std::getline(file, line); // skip LOOKUP_TABLE line

    for (int i = 0; i < nParticles; ++i)
    {
        int fixed_int;
        if (!(file >> fixed_int))
        {
            std::cerr << "Error reading fixed flag for particle " << i << "\n";
            return;
        }
        particles[i].fixed = (fixed_int != 0);
    }

    file.close();
}
