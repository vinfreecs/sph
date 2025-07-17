using sph_float = float;
sph_float r_e = 0.04; // assuming 6 times of the radius
sph_float r_e_b = 0.02;
sph_float PI = 3.14;
sph_float dynamic_viscosity = 10;
sph_float gravity = 0;
sph_float radius = 0.01;
sph_float k = 10;
sph_float rho_o = 998;
sph_float rest_pressure = 0;
int boundary = 2;

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

void boundary_navie(particle &p, double bounce_factor = 0.25)
{
    double hole_size = 0.5;
    double half_boundary = boundary / 2.0;
    double half_hole = hole_size / 2.0;

    // X-axis
    if (p.pos_x - radius < -half_boundary)
    {
        p.pos_x = -half_boundary + radius;
        p.vel_x *= -bounce_factor;
    }

    else if (p.pos_x + radius > half_boundary)
    {
        p.pos_x = half_boundary - radius;
        p.vel_x *= -bounce_factor;
    }

    // Y-axis
    if (p.pos_y - radius < -half_boundary)
    {
        p.pos_y = -half_boundary + radius;
        p.vel_y *= -bounce_factor;
    }
    else if (p.pos_y + radius > half_boundary)
    {
        p.pos_y = half_boundary - radius;
        p.vel_y *= -bounce_factor;
    }
}

sph_float distance_squared(particle p, particle q)
{
    sph_float dx = p.pos_x - q.pos_x;
    sph_float dy = p.pos_y - q.pos_y;
    sph_float r_dist = dx * dx + dy * dy;
    return r_dist;
}

sph_float smoothing_pressure(sph_float r)
{
    if (r < r_e)
    {
        sph_float val = 45 * std::pow(r_e - r, 3) / (PI * std::pow(r_e, 6));
        return val;
    }
    else
        return 0.0;
}

sph_float smoothing_viscosity(sph_float r)
{
    if (r < r_e)
    {
        sph_float val = 45 * (r_e - r) / (PI * std::pow(r_e, 6));
        return val;
    }
    else
        return 0.0;
}

sph_float smoothing_density(sph_float r)
{
    r = std::sqrt(r);
    if (r < r_e)
    {
        sph_float val = 315 * std::pow((r_e * r_e - r * r), 3) / (64 * PI * std::pow(r_e, 9));
        return val;
    }
    else
        return 0.0;
}

void pressure_viscosity_force(particle &p, particle &q, sph_float dt)
{
    sph_float r_squared = distance_squared(p, q);
    sph_float r = std::sqrt(r_squared);
    if (r > r_e)
        return; // no interaction if distance is greater than smoothing length
    // if (r < 1e-6)
    //     r = 1e-6; // to avoid division by zero
    // Calculate the force due to pressure and viscosity

    
    if (q.fixed)
    {
        sph_float min_distance = 2.0f * radius; // assuming both particles have the same radius
        sph_float actual_distance = std::sqrt(distance_squared(p, q));

        if (actual_distance < min_distance && actual_distance > 1e-6f)
        {
            // Normalize direction vector from q to p
            sph_float dx = p.pos_x - q.pos_x;
            sph_float dy = p.pos_y - q.pos_y;
            sph_float inv_dist = 1.0f / actual_distance;
            dx *= inv_dist;
            dy *= inv_dist;

            // Push particle outside the boundary
            sph_float penetration = min_distance - actual_distance;
            p.pos_x += dx * penetration;
            p.pos_y += dy * penetration;

            // Reflect velocity
            sph_float v_dot_n = p.vel_x * dx + p.vel_y * dy;
            p.vel_x -= 2.0f * v_dot_n * dx;
            p.vel_y -= 2.0f * v_dot_n * dy;

            // // Optionally apply damping or restitution
            // float damping = 0.80f; // energy loss on collision
            // p.vel_x *= damping;
            // p.vel_y *= damping;
        }
    }

    if (!q.fixed)
    {
        sph_float w_r_pressure = smoothing_pressure(r);
    sph_float w_r_viscosity = smoothing_viscosity(r);
    sph_float force_pressure_x = q.mass * (p.pressure + q.pressure) * w_r_pressure * (q.pos_x - p.pos_x) / (2 * q.density * r); // TODO
    sph_float force_pressure_y = q.mass * (p.pressure + q.pressure) * w_r_pressure * (q.pos_y - p.pos_y) / (2 * q.density * r);
    sph_float force_viscosity_x = q.mass * (q.vel_x - p.vel_x) * w_r_viscosity / (q.density);
    sph_float force_viscosity_y = q.mass * (q.vel_y - p.vel_y) * w_r_viscosity / (q.density);
    sph_float fx = -force_pressure_x + dynamic_viscosity * force_viscosity_x;
    sph_float fy = -force_pressure_y + dynamic_viscosity * force_viscosity_y;
    p.force_x += fx;
    p.force_y += fy;
    p.vel_x += p.force_x * dt / p.density;
    p.vel_y += p.force_y * dt / p.density;

        q.force_x -= fx;
        q.force_y -= fy;
        q.vel_x += q.force_x * dt / q.density;
        q.vel_y += q.force_y * dt / q.density;
    }
}

inline int get_hash(int x, int y)
{
    std::size_t h = x * 73856093 ^ y * 19349663;
    return h;
}

inline void update_neighbors(std::vector<particle> &particles,
                             std::vector<int> &neighbours,
                             std::vector<int> &particle_neighbours)
{
    int n = particles.size();
    std::fill(neighbours.begin(), neighbours.end(), -1);
    std::fill(particle_neighbours.begin(), particle_neighbours.end(), -1);

    for (int i = 0; i < n; ++i)
    {
        particle &p = particles[i];
        int x = (int)((p.pos_x + (sph_float)boundary / 2.0f) / r_e);
        int y = (int)((p.pos_y + (sph_float)boundary / 2.0f) / r_e);
        int hash = get_hash(x, y);
        int cell_id = ((hash % n) + n) % n; // use a fixed-size table (size = n here)
        int temp = neighbours[cell_id];
        neighbours[cell_id] = i; // atomic_exchange should be used here
        particle_neighbours[i] = temp;
    }
}

void calculation_density_pressure_hash(std::vector<particle> &particles,
                                       const std::vector<int> &neighbours,
                                       const std::vector<int> &particle_neighbours)
{
    int n = particles.size();
    for (int i = 0; i < n; i++)
    {
        particle &p = particles[i];
        p.density = 999.0f;
        if (!p.fixed)
        {
            p.density = 0;
            p.pressure = 0;
            int x = (int)((p.pos_x + (sph_float)boundary / 2.0f) / r_e);
            int y = (int)((p.pos_y + (sph_float)boundary / 2.0f) / r_e);
            for (int dx = -1; dx <= 1; ++dx)
            {
                for (int dy = -1; dy <= 1; ++dy)
                {
                    int hx = x + dx;
                    int hy = y + dy;
                    int hash = get_hash(hx, hy);
                    int cell_id = ((hash % n) + n) % n; // same fix as above

                    int j = neighbours[cell_id];
                    while (j != -1)
                    {
                        if (j != i || j == i)
                        {
                            particle &q = particles[j];
                            sph_float r = distance_squared(p, q);
                            sph_float w_r = smoothing_density(r);
                            p.density += q.mass * (w_r);
                        }
                        j = particle_neighbours[j];
                    }
                }
            }
        }
        // if (p.density < 1e-5) p.density = 1e-5;
        p.pressure = rest_pressure + k * (p.density - rho_o);
    }
}

void calculate_force_hash(std::vector<particle> &particles,
                          std::vector<int> &neighbours,
                          std::vector<int> &particle_neighbours, sph_float dt)
{
    int n = particles.size();
    for (int i = 0; i < n; i++)
    {
        particle &p = particles[i];
        if (!p.fixed)
        {

            int x = (int)((p.pos_x + (sph_float)boundary / 2.0f) / r_e);
            int y = (int)((p.pos_y + (sph_float)boundary / 2.0f) / r_e);
            for (int dx = -1; dx <= 1; ++dx)
            {
                for (int dy = -1; dy <= 1; ++dy)
                {
                    int hx = x + dx;
                    int hy = y + dy;
                    int hash = get_hash(hx, hy);
                    int cell_id = ((hash % n) + n) % n; // same fix as above

                    int j = neighbours[cell_id];
                    while (j != -1)
                    {
                        if (i < j)
                        {
                            particle &q = particles[j];

                            sph_float r = distance_squared(p, q);
                            if (r < r_e)

                            {
                                pressure_viscosity_force(p, q, dt);
                            }
                        }
                        j = particle_neighbours[j];
                    }
                }
            }
        }
        p.force_y -= p.density * gravity; // gravity only in the -ve z direction
        p.vel_y += p.force_y * dt / p.density;
    }
}

void reset_forces(std::vector<particle> &particles)
{
    for (auto &p : particles)
    {
        if (!p.fixed)
        {
            p.force_x = 0;
            p.force_y = 0;
        }
    }
}

void update_positions(std::vector<particle> &particles, sph_float dt)
{
    for (auto &p : particles)
    {
        if (!p.fixed)
        {
            p.pos_x += p.vel_x * dt;
            p.pos_y += p.vel_y * dt;
            boundary_navie(p);
        }
    }
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
        file >> particles[i].pos_x >> particles[i].pos_y;
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
        file >> particles[i].vel_x >> particles[i].vel_y;
    }

    for (int i = 0; i < nParticles; i++)
    {
        particles[i].id = i;
        particles[i].radius = radius;
    }

    file.close();
}

void write_particles_square_format(
    const std::string &filename,
    int total_particles,
    double square_size,
    double mass = 1.0,
    double vx = 0.0, double vy = 0.0)
{
    std::ofstream file(filename);
    if (!file.is_open())
    {
        std::cerr << "Error: Could not open file " << filename << "\n";
        return;
    }

    int particles_per_axis = static_cast<int>(std::ceil(std::sqrt(total_particles)));
    double spacing = square_size / particles_per_axis;

    // Compute offset to center vertically (Y-axis centered at 0)
    double y_offset = -square_size / 2.0;

    // Write positions
    file << "positions\n";
    int particles_written = 0;
    for (int i = 0; i < particles_per_axis && particles_written < total_particles; ++i)
    {
        for (int j = 0; j < particles_per_axis && particles_written < total_particles; ++j)
        {
            double x = i * spacing + spacing / 2.0 + y_offset;
            double y = j * spacing + spacing / 2.0 + y_offset; // Center Y
            file << x << " " << y << "\n";
            particles_written++;
        }
    }

    // Write masses
    file << "masses\n";
    for (int i = 0; i < total_particles; ++i)
    {
        file << mass << "\n";
    }

    // Write velocities
    file << "velocities\n";
    for (int i = 0; i < total_particles; ++i)
    {
        file << vx << " " << vy << "\n";
    }

    file.close();
    std::cout << "Wrote " << total_particles << " particles to " << filename << "\n";
}

std::vector<particle> create_boundary_particles_2d(sph_float boundary_size, sph_float r_e_b, float overlap_factor = 0.2f)
{
    std::vector<particle> boundary_particles;

    sph_float spacing = r_e_b * overlap_factor; // reduce spacing to overlap
    int num_per_axis = boundary_size / spacing;
    sph_float half_size = boundary_size / 2.0f;

    for (int i = 0; i <= num_per_axis; ++i)
    {
        for (int j = 0; j <= num_per_axis; ++j)
        {
            sph_float x = -half_size + i * spacing;
            sph_float y = -half_size + j * spacing;

            bool is_border = (i == 0 || i == num_per_axis || j == 0 || j == num_per_axis);

            if (is_border)
            {
                particle p;
                p.pos_x = x;
                p.pos_y = y;
                p.vel_x = 0.0f;
                p.vel_y = 0.0f;
                p.mass = 1.0f;
                p.density = 7870.0f; // optional: high to simulate solid wall
                p.pressure = 0.0f;
                p.force_x = 0.0f;
                p.force_y = 0.0f;
                p.fixed = true;
                p.radius = r_e_b / 2.0f;
                boundary_particles.push_back(p);
            }
        }
    }

    return boundary_particles;
}

std::vector<particle> create_boundary_circle_2d(sph_float radius, sph_float r_e_b, float overlap_factor = 0.2f)
{
    std::vector<particle> boundary_particles;

    sph_float spacing = r_e_b * overlap_factor;
    int num_particles = static_cast<int>((2.0f * M_PI * radius) / spacing); // circumference / spacing
    sph_float angle_step = 2.0f * M_PI / num_particles;

    for (int i = 0; i < num_particles; ++i)
    {
        sph_float angle = i * angle_step;
        sph_float x = radius * std::cos(angle);
        sph_float y = radius * std::sin(angle);

        particle p;
        p.pos_x = x;
        p.pos_y = y;
        p.vel_x = 0.0f;
        p.vel_y = 0.0f;
        p.mass = 1.0f;
        p.density = 7870.0f; // optional: solid wall
        p.pressure = 0.0f;
        p.force_x = 0.0f;
        p.force_y = 0.0f;
        p.fixed = true;
        p.radius = r_e_b / 2.0f;
        boundary_particles.push_back(p);
    }

    return boundary_particles;
}

std::vector<particle> create_boundary_circle_2d_with_opening(
    sph_float radius,
    sph_float r_e_b,
    float overlap_factor = 0.2f,
    sph_float opening_angle_rad = M_PI / 6.0, // 30° gap
    sph_float opening_center_rad = 0.0        // gap centered at angle 0 (right side)
)
{
    std::vector<particle> boundary_particles;

    sph_float spacing = r_e_b * overlap_factor;
    int num_particles = static_cast<int>((2.0f * M_PI * radius) / spacing);
    sph_float angle_step = 2.0f * M_PI / num_particles;

    sph_float half_opening = opening_angle_rad / 2.0;

    for (int i = 0; i < num_particles; ++i)
    {
        sph_float angle = i * angle_step;

        // Check if angle is within opening range (gap)
        sph_float delta_angle = std::fmod(angle - opening_center_rad + 2 * M_PI, 2 * M_PI);
        if (delta_angle < opening_angle_rad)
        {
            continue; // skip this particle to create the gap
        }

        sph_float x = radius * std::cos(angle);
        sph_float y = radius * std::sin(angle);

        particle p;
        p.pos_x = x;
        p.pos_y = y;
        p.vel_x = 0.0f;
        p.vel_y = 0.0f;
        p.mass = 1.0f;
        p.density = 7870.0f;
        p.pressure = 0.0f;
        p.force_x = 0.0f;
        p.force_y = 0.0f;
        p.fixed = true;
        p.radius = r_e_b / 2.0f;

        boundary_particles.push_back(p);
    }

    return boundary_particles;
}