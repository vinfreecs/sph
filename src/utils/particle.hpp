class particle{
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
                double x = i * spacing + spacing/2;  // Center particles in their cells
                double y = j * spacing + spacing/2;
                double z = k * spacing + spacing/2;
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

std::vector<particle> create_boundary_particles(sph_float boundary_size, sph_float r_e, sph_float wall_thickness = 1.0) {
    std::vector<particle> boundary_particles;

    int num_per_axis = boundary_size / r_e;
    int thickness_layers = wall_thickness / r_e;

    for (int i = 0; i <= num_per_axis; ++i) {
        for (int j = 0; j <= num_per_axis; ++j) {
            for (int k = 0; k <= num_per_axis; ++k) {
                sph_float x = -boundary_size / 2 + i * r_e;
                sph_float y = -boundary_size / 2 + j * r_e;
                sph_float z = -boundary_size / 2 + k * r_e;

                bool is_on_wall = 
                    (i < thickness_layers) || (i > num_per_axis - thickness_layers) ||
                    (j < thickness_layers) || (j > num_per_axis - thickness_layers) ||
                    (k < thickness_layers) || (k > num_per_axis - thickness_layers);

                if (is_on_wall) {
                    particle p;
                    p.pos_x = x;
                    p.pos_y = y;
                    p.pos_z = z;
                    p.vel_x = 0.0f;
                    p.vel_y = 0.0f;
                    p.vel_z = 0.0f;
                    p.mass = 1.0f;  // Can be tuned
                    p.density = 0.0f;
                    p.pressure = 0.0f;
                    p.force_x = p.force_y = p.force_z = 0.0f;
                    p.fixed = true;
                    boundary_particles.push_back(p);
                }
            }
        }
    }

    return boundary_particles;
}
