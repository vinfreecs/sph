class particle{
public:
    int id;
    int bucket_id=-1; // For neighbour search
    float pos_x, pos_y, pos_z;
    float vel_x, vel_y, vel_z;
    float mass;
    float density;
    float pressure;
    float force_x, force_y, force_z;
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

