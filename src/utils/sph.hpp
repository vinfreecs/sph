sph_float distance_squared(particle p, particle q)
{
    sph_float dx = p.pos_x - q.pos_x;
    sph_float dy = p.pos_y - q.pos_y;
    sph_float dz = p.pos_z - q.pos_z;
    sph_float r_dist = dx * dx + dy * dy + dz * dz;
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

    // Z-axis
    if (p.pos_z - radius < -half_boundary)
    {
        p.pos_z = -half_boundary + radius;
        p.vel_z *= -bounce_factor;
    }
    else if (p.pos_z + radius > half_boundary)
    {
        p.pos_z = half_boundary - radius;
        p.vel_z *= -bounce_factor;
    }
}

void pressure_viscosity_force(particle &p, particle &q)
{
    sph_float r_squared = distance_squared(p, q);
    sph_float r = std::sqrt(r_squared);
    if (r > r_e)
        return; // no interaction if distance is greater than smoothing length
    // if (r < 1e-5)
    //     r = 1e-5; // to avoid division by zero
    // Calculate the force due to pressure and viscosity
    if (q.fixed)
    {
        sph_float min_distance = p.radius + q.radius + 0.01; // assuming both particles have the same radius
        sph_float actual_distance = std::sqrt(distance_squared(p, q));

        if (actual_distance <= min_distance)
        {
            // Normalize direction vector from q to p
            sph_float dx = p.pos_x - q.pos_x;
            sph_float dy = p.pos_y - q.pos_y;
            sph_float dz = p.pos_z - q.pos_z;

            sph_float inv_dist = 1.0f / actual_distance;
            dx *= inv_dist;
            dy *= inv_dist;
            dz *= inv_dist;

            // Push particle outside the boundary
            sph_float penetration = min_distance - actual_distance;
            p.pos_x += dx * penetration;
            p.pos_y += dy * penetration;
            p.pos_z += dz * penetration;

            // Reflect velocity
            sph_float v_dot_n = p.vel_x * dx + p.vel_y * dy + p.vel_z * dz;
            p.vel_x -= 1.0f * v_dot_n * dx;
            p.vel_y -= 1.0f * v_dot_n * dy;
            p.vel_z -= 1.0f * v_dot_n * dz;

            // // Optionally apply damping or restitution
            // float damping = 1.0f; // energy loss on collision
            // p.vel_x *= damping;
            // p.vel_y *= damping;
        }
    }
    sph_float w_r_pressure = smoothing_pressure(r);
    sph_float w_r_viscosity = smoothing_viscosity(r);
    sph_float force_pressure_x = q.mass * (p.pressure + q.pressure) * w_r_pressure * (q.pos_x - p.pos_x) / (2 * q.density * r); // TODO
    sph_float force_pressure_y = q.mass * (p.pressure + q.pressure) * w_r_pressure * (q.pos_y - p.pos_y) / (2 * q.density * r);
    sph_float force_pressure_z = q.mass * (p.pressure + q.pressure) * w_r_pressure * (q.pos_z - p.pos_z) / (2 * q.density * r);
    sph_float force_viscosity_x = q.mass * (q.vel_x - p.vel_x) * w_r_viscosity / (q.density);
    sph_float force_viscosity_y = q.mass * (q.vel_y - p.vel_y) * w_r_viscosity / (q.density);
    sph_float force_viscosity_z = q.mass * (q.vel_z - p.vel_z) * w_r_viscosity / (q.density);
    sph_float fx = -force_pressure_x + dynamic_viscosity * force_viscosity_x;
    sph_float fy = -force_pressure_y + dynamic_viscosity * force_viscosity_y;
    sph_float fz = -force_pressure_z + dynamic_viscosity * force_viscosity_z;
    p.force_x += fx;
    p.force_y += fy;
    p.force_z += fz;
    if (!q.fixed)
    {
        q.force_x -= fx;
        q.force_y -= fy;
        q.force_z -= fz;
    }
}

inline int get_hash(int cell_x, int cell_y, int cell_z)
{
    return (cell_x * 9283711) ^ (cell_y * 689287499) ^ (cell_z * 283923481);
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
        int z = (int)((p.pos_z + (sph_float)boundary / 2.0f) / r_e);
        int hash = get_hash(x, y, z);
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
            int z = (int)((p.pos_z + (sph_float)boundary / 2.0f) / r_e);
            for (int dx = -1; dx <= 1; ++dx)
            {
                for (int dy = -1; dy <= 1; ++dy)
                {
                    for (int dz = -1; dz <= 1; ++dz)
                    {
                        int hx = x + dx;
                        int hy = y + dy;
                        int hz = z + dz;
                        int hash = get_hash(hx, hy, hz);
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
        }
        // if (p.density < 1e-5) p.density = 1e-5;
        p.pressure = rest_pressure + k * (p.density - rho_o);
    }
}

void calculate_force_hash(std::vector<particle> &particles,
                          std::vector<int> &neighbours,
                          std::vector<int> &particle_neighbours)
{
    int n = particles.size();
    for (int i = 0; i < n; i++)
    {
        particle &p = particles[i];
        if (!p.fixed)
        {

            int x = (int)((p.pos_x + (sph_float)boundary / 2.0f) / r_e);
            int y = (int)((p.pos_y + (sph_float)boundary / 2.0f) / r_e);
            int z = (int)((p.pos_z + (sph_float)boundary / 2.0f) / r_e);
            for (int dx = -1; dx <= 1; ++dx)
            {
                for (int dy = -1; dy <= 1; ++dy)
                {
                    for (int dz = -1; dz <= 1; ++dz)
                    {
                        int hx = x + dx;
                        int hy = y + dy;
                        int hz = z + dz;
                        int hash = get_hash(hx, hy, hz);
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
                                    pressure_viscosity_force(p, q);
                                }
                            }
                            j = particle_neighbours[j];
                        }
                    }
                }
            }
        }
        p.force_y -= p.density * gravity; // gravity only in the -ve z direction
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
            p.force_z = 0;
        }
    }
}

void update_velocity(std::vector<particle> &particles, sph_float dt)
{
    for (auto &p : particles)
    {
        // std::cout <<"Force x :"<<p.force_x<<" y :"<<p.force_y<<" z :"<<p.force_z<<std::endl;
        if (!p.fixed)
        {

            p.vel_x += p.force_x * dt / p.density;
            p.vel_y += p.force_y * dt / p.density;
            p.vel_z += ((p.force_z / p.density)) * dt;
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
            p.pos_z += p.vel_z * dt;
            boundary_navie(p);
        }
    }
}
