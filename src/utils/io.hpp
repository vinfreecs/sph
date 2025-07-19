std::string generateVTKFilename(int frame)
{
    return "frame_" + std::to_string(frame) + ".vtk";
}

// void writeVTK(const std::string &filename, std::vector<particle> &particles)
// {
//     std::ofstream out(filename);
//     if (!out.is_open())
//     {
//         std::cerr << "Failed to open file: " << filename << std::endl;
//         return;
//     }

//     size_t numParticles = 0;
//     for (auto &p : particles)
//         if (!p.fixed)
//             ++numParticles;

//     out << "# vtk DataFile Version 4.0\n";
//     out << "hesp visualization file\n";
//     out << "ASCII\n";
//     out << "DATASET UNSTRUCTURED_GRID\n";
//     out << "POINTS " << numParticles << " double\n";

//     out << std::fixed << std::setprecision(6);
//     for (auto particle : particles)
//     {
//         if (particle.fixed)
//             continue;
//         double x = particle.pos_x;
//         double y = particle.pos_y;
//         double z = particle.pos_z;
//         out << x << " " << y << " " << z << "\n";
//     }

//     out << "CELLS 0 0\n";
//     out << "CELL_TYPES 0\n";
//     out << "POINT_DATA " << numParticles << "\n";

//     // Write masses
//     out << "SCALARS m double\n";
//     out << "LOOKUP_TABLE default\n";
//     for (auto m : particles)
//     {
//         if (m.fixed)
//             continue;
//         out << m.mass << "\n";
//     }

//     // Write velocities
//     out << "VECTORS v double\n";
//     for (auto particle : particles)
//     {
//         if (particle.fixed)
//             continue;
//         double vx = particle.vel_x;
//         double vy = particle.vel_y;
//         double vz = particle.vel_z;
//         out << vx << " " << vy << " " << vz << "\n";
//     }

//     out.close();
// }

void writeVTK(const std::string &filename, std::vector<particle> &particles)
{
    std::ofstream out(filename);
    if (!out.is_open())
    {
        std::cerr << "Failed to open file: " << filename << std::endl;
        return;
    }

    size_t numParticles = 0;
    for (const auto &p : particles)
        if (!p.fixed)
            ++numParticles;

    out << "# vtk DataFile Version 4.0\n";
    out << "SPH particle data\n";
    out << "ASCII\n";
    out << "DATASET UNSTRUCTURED_GRID\n";
    out << "POINTS " << numParticles << " double\n";

    out << std::fixed << std::setprecision(6);
    for (const auto &p : particles)
    {
        if (p.fixed) continue;
        out << p.pos_x << " " << p.pos_y << " " << p.pos_z << "\n";
    }

    out << "CELLS 0 0\n";
    out << "CELL_TYPES 0\n";
    out << "POINT_DATA " << numParticles << "\n";

    // Scalar: radius
    out << "SCALARS radius double\n";
    out << "LOOKUP_TABLE default\n";
    for (const auto &p : particles)
        if (!p.fixed)
            out << p.radius << "\n";

    // Scalar: mass
    out << "SCALARS mass double\n";
    out << "LOOKUP_TABLE default\n";
    for (const auto &p : particles)
        if (!p.fixed)
            out << p.mass << "\n";

    // Scalar: density
    out << "SCALARS density double\n";
    out << "LOOKUP_TABLE default\n";
    for (const auto &p : particles)
        if (!p.fixed)
            out << p.density << "\n";

    // Scalar: pressure
    out << "SCALARS pressure double\n";
    out << "LOOKUP_TABLE default\n";
    for (const auto &p : particles)
        if (!p.fixed)
            out << p.pressure << "\n";

    // Vector: velocity
    out << "VECTORS velocity double\n";
    for (const auto &p : particles)
        if (!p.fixed)
            out << p.vel_x << " " << p.vel_y << " " << p.vel_z << "\n";

    // Vector: force
    out << "VECTORS force double\n";
    for (const auto &p : particles)
        if (!p.fixed)
            out << p.force_x << " " << p.force_y << " " << p.force_z << "\n";

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

    for (int i = 0; i < nParticles; i++)
    {
        particles[i].id = i;
        particles[i].radius = radius_device;
        particles[i].fixed = false;
    }

    file.close();
}
