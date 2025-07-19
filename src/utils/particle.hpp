class particle
{
public:
    int id;
    sph_float radius;
    sph_float pos_x, pos_y, pos_z;
    sph_float vel_x, vel_y, vel_z;
    sph_float mass;
    sph_float density;
    sph_float pressure;
    sph_float force_x, force_y, force_z;
    bool fixed;
};

void write_particles_cube_format(
    const std::string& filename,
    int total_particles,
    double cube_size,
    double mass = 1.0,
    double vx = 0.0, double vy = 0.0, double vz = 0.0
) {
    std::ofstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error: Could not open file " << filename << "\n";
        return;
    }

    // Calculate how many particles fit along each axis
    int particles_per_axis = static_cast<int>(std::ceil(std::cbrt(total_particles)));
    double spacing = cube_size / particles_per_axis;

    // Write positions
    file << "positions\n";
    int particles_written = 0;
    for (int i = 0; i < particles_per_axis && particles_written < total_particles; ++i) {
        for (int j = 0; j < particles_per_axis && particles_written < total_particles; ++j) {
            for (int k = 0; k < particles_per_axis && particles_written < total_particles; ++k) {
                double x = i * spacing + spacing/2 - cube_size/2;  // Center particles in their cells
                double y = j * spacing + spacing/2 - cube_size/2;
                double z = k * spacing + spacing/2 + cube_size/2;
                file << x << " " << y << " " << z << "\n";
                particles_written++;
            }
        }
    }

    // Write masses
    file << "masses\n";
    for (int i = 0; i < total_particles; ++i) {
        file << mass << "\n";
    }

    // Write velocities
    file << "velocities\n";
    for (int i = 0; i < total_particles; ++i) {
        file << vx << " " << vy << " " << vz << "\n";
    }

    file.close();
    std::cout << "Wrote " << total_particles << " particles to " << filename << "\n";
}


std::vector<particle> create_boundary_particles_3d(sph_float boundary_size, sph_float r_e_b, float overlap_factor = 0.2f)
{
    std::vector<particle> boundary_particles;

    sph_float spacing = r_e_b * overlap_factor; // reduce spacing to overlap
    int num_per_axis = boundary_size / spacing;
    sph_float half_size = boundary_size / 2.0f;

    for (int i = 0; i <= num_per_axis; ++i)
    {
        for (int j = 0; j <= num_per_axis; ++j)
        {
            for (int k = 0; k <= num_per_axis; ++k)
            {
                sph_float x = -half_size + i * spacing;
                sph_float y = -half_size + j * spacing;
                sph_float z = -half_size + k * spacing;

                // If the particle is on any boundary face
                bool is_border = (
                    i == 0 || i == num_per_axis ||
                    j == 0 || j == num_per_axis ||
                    k == 0 || k == num_per_axis
                );

                if (is_border)
                {
                    particle p;
                    p.pos_x = x;
                    p.pos_y = y;
                    p.pos_z = z;
                    p.vel_x = 0.0f;
                    p.vel_y = 0.0f;
                    p.vel_z = 0.0f;
                    p.mass = 1.0f;
                    p.density = 7870.0f;
                    p.pressure = 0.0f;
                    p.force_x = 0.0f;
                    p.force_y = 0.0f;
                    p.force_z = 0.0f;
                    p.fixed = true;
                    p.radius = r_e_b / 2.0f;
                    boundary_particles.push_back(p);
                }
            }
        }
    }

    return boundary_particles;
}

std::vector<particle> create_fluid_sphere_particles_3d(
    sph_float sphere_radius, sph_float r_e)
{
    std::vector<particle> fluid_particles;

    sph_float spacing = r_e;  // Proper spacing for fluid
    int num_per_axis = static_cast<int>((2.0f * sphere_radius) / spacing);
    sph_float half_extent = num_per_axis * spacing / 2.0f;

    for (int i = 0; i <= num_per_axis; ++i)
    {
        for (int j = 0; j <= num_per_axis; ++j)
        {
            for (int k = 0; k <= num_per_axis; ++k)
            {
                sph_float x = -half_extent + i * spacing;
                sph_float y = -half_extent + j * spacing;
                sph_float z = -half_extent + k * spacing;

                // Only keep particles within the sphere radius
                if (x*x + y*y + z*z <= sphere_radius * sphere_radius)
                {
                    particle p;
                    p.pos_x = x;
                    p.pos_y = y;
                    p.pos_z = z;
                    p.vel_x = 0.0f;
                    p.vel_y = 0.0f;
                    p.vel_z = 0.0f;
                    p.mass = 1.0f;
                    p.density = 998.29f; // typical density of water
                    p.pressure = 0.0f;
                    p.force_x = 0.0f;
                    p.force_y = 0.0f;
                    p.force_z = 0.0f;
                    p.fixed = false; // it's a fluid particle
                    p.radius = r_e / 2.0f; // or whatever makes sense in your kernel
                    fluid_particles.push_back(p);
                }
            }
        }
    }

    return fluid_particles;
}
