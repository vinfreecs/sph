#include <vector>
#include <iostream>

using sph_float = float;
float PI = 3.14;
float r_e = 0.05;
sph_float dynamic_viscosity = 0.001;

struct particle
{
    int id;
    float x, y;
    float radius;
    float mass;
    float pressure;
    float density;
    float velocity_x, velocity_y;
    int bucket_id = -1;
    bool is_boundary = false;
    // particle type fluid and boundary
};



std::vector<particle> create_boundary_particles(
    float box_width,
    float box_height,
    float spacing,
    float radius,
    float mass,
    int layers = 1)  // Number of ghost layers
{
    std::vector<particle> boundaries;
    int id_counter = 0;

    int fluid_cols = static_cast<int>(box_width / spacing);
    int fluid_rows = static_cast<int>(box_height / spacing);

    float offset_x = -0.5f * box_width;
    float offset_y = -0.5f * box_height;

    int total_cols = fluid_cols + 2 * layers;
    int total_rows = fluid_rows + 2 * layers;

    for (int i = 0; i < total_rows; ++i) {
        for (int j = 0; j < total_cols; ++j) {
            // Inside fluid region → skip
            if (i >= layers && i < (total_rows - layers) &&
                j >= layers && j < (total_cols - layers)) {
                continue;
            }

            particle p;
            p.x = offset_x - layers * spacing + j * spacing;
            p.y = offset_y - layers * spacing + i * spacing;
            p.radius = radius;
            p.mass = mass;
            p.density = 998.0f;
            p.pressure = 0.0f;
            p.velocity_x = 0.0f;
            p.velocity_y = 0.0f;
            p.bucket_id = -1;
            p.is_boundary = true;
            p.id = id_counter++;
            boundaries.push_back(p);
        }
    }

    return boundaries;
}



std::vector<particle> create_particles_grid(int rows, int cols, float spacing, float radius, float mass)
{
    std::vector<particle> particles;
    float offset_x = -0.5f * (cols - 1) * spacing;
    float offset_y = -0.5f * (rows - 1) * spacing;

    for (int i = 0; i < rows; ++i)
    {
        for (int j = 0; j < cols; ++j)
        {
            particle p;
            p.x = offset_x + j * spacing;
            p.y = offset_y + i * spacing;
            p.radius = radius;
            p.mass = mass;
            p.density = 0.0f;
            p.pressure = 0.0f;
            p.velocity_x = 0.0f;
            p.velocity_y = 0.0f;
            p.id = i * rows + j;
            p.bucket_id = -1;
            particles.push_back(p);
        }
    }
    return particles;
}

std::vector<particle> create_all_particles()
{
    std::vector<particle> fluid_particles = create_particles_grid(
        10, 10,
        0.04f,
        0.01f,
        1.0f);

    std::vector<particle> boundary_particles = create_boundary_particles(
        1.5f, 1.5f, // width, height of the container box
        0.001f,     // spacing
        0.002f,     // radius
        1.0f        // mass
    );
    std::vector<particle> particles;
    particles.insert(particles.end(), fluid_particles.begin(), fluid_particles.end());
    particles.insert(particles.end(), boundary_particles.begin(), boundary_particles.end());

    return particles;
}

// Combine them
std::vector<particle> particles = create_all_particles();

bool is_near_box_boundary(const particle& p, float box_half_size = 1.0f, float threshold = 0.05f) {
    float left   = -box_half_size;
    float right  =  box_half_size;
    float bottom = -box_half_size;
    float top    =  box_half_size;

    return (
        std::abs(p.x - left)   < threshold ||
        std::abs(p.x - right)  < threshold ||
        std::abs(p.y - bottom) < threshold ||
        std::abs(p.y - top)    < threshold
    );
}

void update_buckets(std::vector<particle> &particles, std::vector<std::vector<int>> &buckets, int boundary, int nBuckets)
{
    for (int i = 0; i < particles.size(); i++)
    {
        sph_float x = particles[i].x + (sph_float)boundary / 2;
        sph_float y = particles[i].y + (sph_float)boundary / 2;
        if (x < 0 || x > boundary || y < 0 || y > boundary)
        {
            std::cerr << "Particle " << i << " is out of bounds. Please check the input file." << std::endl;
            std::cerr << "x: " << x << " y: " << y << " z: " << " bucket_idx: " << std::endl;
            std::cerr << "nBuckets: " << nBuckets << " boundary: " << boundary << std::endl;
            continue;
        }
        int bucket_idx = (int)(x / r_e) + int(y / r_e) * nBuckets;

        if (bucket_idx > nBuckets * nBuckets * nBuckets - 1 || bucket_idx < 0)
        {
            std::cerr << "Particle " << i << " is out of bounds. Please check the input file." << std::endl;
            std::cerr << "x: " << x << " y: " << y << " z: " << " bucket_idx: " << bucket_idx << std::endl;
            std::cerr << "nBuckets: " << nBuckets << " boundary: " << boundary << std::endl;
            continue;
        }
        if (particles[i].bucket_id == bucket_idx)
        {
            continue; // already in the same bucket
        }
        else
        {
            if (particles[i].bucket_id != -1)
            {
                // Removing from old bucket
                int old_bucket_idx = particles[i].bucket_id;
                for (int j = 0; j < buckets[old_bucket_idx].size(); j++)
                {
                    if (buckets[old_bucket_idx][j] == particles[i].id)
                    {
                        buckets[old_bucket_idx].erase(buckets[old_bucket_idx].begin() + j);
                        break;
                    }
                }
            }
            particles[i].bucket_id = bucket_idx;
            buckets[bucket_idx].push_back(particles[i].id);
        }

        // std::cout<<"Particle "<<i<<" in bucket "<<bucket_idx<<std::endl;
    }
}

int boundary = 3;
int nBuckets = ceil(boundary / r_e);

std::vector<std::vector<int>> buckets(nBuckets *nBuckets); // neighbours 3 is the boundary and 0.006 is the smoothing length effective radius etc

std::vector<int> neighbour_list(particles.size());
std::vector<int> start_neighbour(particles.size(), -1);

int calculate_cell_hash(particle &p)
{
    int x = std::floor(p.x / r_e);
    int y = std::floor(p.x / r_e);

    int x_hash = x * 73856093;
    int y_hash = y * 19349663;
    return x_hash + y_hash;
    // 83,492,791 for z
}

void update_neighbours_list(std::vector<int> neighbour_list, std::vector<particle> particles)
{
    for (int i = 0; i < particles.size(); i++)
    {
        particle &p = particles[i];
        int cell_ind = calculate_cell_hash(p) % (neighbour_list.size());
        // neighbour_list[i] =
    }
}

void apply_gravity(std::vector<particle> &particles, float gravity, float dt)
{
    for (auto &p : particles)
    {
        if (!p.is_boundary)
            p.velocity_y += gravity * dt; // gravity affects only y-direction
    }
}

void update_positions(std::vector<particle> &particles, float dt)
{
    for (auto &p : particles)
    {
        if (!p.is_boundary)
        {
            p.x += p.velocity_x * dt;
            p.y += p.velocity_y * dt;

    std::cout<<"Particles x "<<p.x<<" y "<<p.y<<std::endl;
    std::cout<<"Particles vel_x "<<p.velocity_x<<" y "<<p.velocity_y<<std::endl;


            
        }
    }
}

void apply_boundary_conditions(std::vector<particle> &particles, float window_width, float window_height, float bounce = -0.25f)
{
    for (auto &p : particles)
    {
        // Left or right wall
        if (p.x - p.radius < -window_width / 2.0f)
        {
            p.x = -window_width / 2.0f + p.radius;
            p.velocity_x *= bounce;
        }
        else if (p.x + p.radius > window_width / 2.0f)
        {
            p.x = window_width / 2.0f - p.radius;
            p.velocity_x *= bounce;
        }

        // Bottom or top wall
        if (p.y - p.radius < -window_height / 2.0f)
        {
            p.y = -window_height / 2.0f + p.radius;
            p.velocity_y *= bounce;
        }
        else if (p.y + p.radius > window_height / 2.0f)
        {
            p.y = window_height / 2.0f - p.radius;
            p.velocity_y *= bounce;
        }
    }
}

float smoothing_density(sph_float r)
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

float distance_squared(particle p, particle q)
{
    sph_float dx = p.x - q.x;
    sph_float dy = p.y - q.y;
    sph_float r_dist = dx * dx + dy * dy;
    return r_dist;
}

void calculate_density(std::vector<particle> &particles)
{
    for (auto &p : particles)
    {
        p.density = 0.0;
        for (auto &q : particles)
        {
            float r_dist_squared = distance_squared(p, q);
            float w_r = smoothing_density(r_dist_squared);
            p.density += q.mass * (w_r);
        }
    }
}

void calculate_pressure(std::vector<particle> &particles, sph_float k = 3, sph_float rho_o = 998.8, sph_float rest_pressure = 10.0)
{
    for (auto &p : particles)
    {
        if (!p.is_boundary)
            p.pressure = rest_pressure + k * (p.density - rho_o);
        // std::cout<<"Particle "<<p.id<<" density :"<<p.density<<" pressure :"<<rest_pressure + k * (p.density - rho_o)<<std::endl;
    }
}

sph_float smoothing_pressure(sph_float r)
{
    if (r < 1e-5)
        r = 1e-5; // to avoid division by zero
    if (r < r_e)
    {
        sph_float val = 45 * std::pow(r_e - r, 3) / (PI * std::pow(r_e, 6) * r);
        return val;
    }
    else
        return 0.0;
}

sph_float smoothing_viscosity(sph_float r)
{
    if (r_e > r)
    {
        sph_float val = 45 * (r_e - r) / (PI * std::pow(r_e, 6));
        return val;
    }
    else
        return 0.0;
}

void pressure_viscosity_force(particle &p, particle &q, float dt)
{
    sph_float r_squared = distance_squared(p, q);
    sph_float r = std::sqrt(r_squared);
    if (r > r_e)
        return; // no interaction if distance is greater than smoothing length
    // if(r < 1e-5) return; // to avoid division by zero
    //  Calculate the force due to pressure and viscosity

    sph_float w_r_pressure = smoothing_pressure(r);
    sph_float w_r_viscosity = smoothing_viscosity(r);
    sph_float force_pressure_x = q.mass * (p.pressure + q.pressure) * w_r_pressure * (q.x - p.x) / (2 * q.density); // TODO
    sph_float force_pressure_y = q.mass * (p.pressure + q.pressure) * w_r_pressure * (q.y - p.y) / (2 * q.density);
    sph_float force_viscosity_x = q.mass * (q.velocity_x - p.velocity_x) * w_r_viscosity / (q.density);
    sph_float force_viscosity_y = q.mass * (q.velocity_y - p.velocity_y) * w_r_viscosity / (q.density);
    sph_float fx = -force_pressure_x;
    sph_float fy = -force_pressure_y;

    p.velocity_x += ((-force_pressure_x / p.density) + (dynamic_viscosity * force_viscosity_x)) * dt;
    p.velocity_y += ((-force_pressure_y / p.density) + (dynamic_viscosity * force_viscosity_y)) * dt;

    if (!q.is_boundary)
    {
        q.velocity_x -= ((-force_pressure_x / p.density) + (dynamic_viscosity * force_viscosity_x)) * dt;
        q.velocity_y -= ((-force_pressure_y / p.density) + (dynamic_viscosity * force_viscosity_y)) * dt;
    }
}

void calculate_force(std::vector<particle> &particles, float dt)
{
    // Force = gravity_force + pressure_force + viscosity_force
    for (int i = 0; i < particles.size(); i++)
    {
        particle &p = particles[i];
        for (int j = 0; j < particles.size(); j++)
        {
            particle &q = particles[j];
            if (i != j)
            {
                pressure_viscosity_force(p, q, dt);
            }
        }
    }
}

std::vector<int> get_neighbours(int idx, int n)
{
    int x = idx % n;
    int y = idx / n;

    std::vector<int> neighbors;
    for (int dx = -1; dx <= 1; dx++)
    {
        for (int dy = -1; dy <= 1; dy++)
        {
            int nx = x + dx;
            int ny = y + dy;
            if (nx >= 0 && ny >= 0 && nx < n && ny < n)
            {
                int neighbor_idx = nx + ny * n;
                neighbors.push_back(neighbor_idx);
            }
        }
    }
    return neighbors;
}

void calculate_density_buckets(std::vector<particle> &particles, std::vector<std::vector<int>> &buckets, int nBuckets)
{
    for (auto &p : particles)
    {
        if (p.is_boundary)
            continue;
        p.density = 0.0;
        int p_bucket_id = p.bucket_id;
        std::vector<int> neighbour_buckets = get_neighbours(p_bucket_id, nBuckets);
        for (int j = 0; j < neighbour_buckets.size(); j++)
        {
            int neighbour_bucket_id = neighbour_buckets[j];
            if (neighbour_bucket_id < 0 || neighbour_bucket_id >= nBuckets * nBuckets * nBuckets)
                continue;
            for (int k = 0; k < buckets[neighbour_bucket_id].size(); k++)
            {
                particle &q = particles[buckets[neighbour_bucket_id][k]];
                sph_float r_dist_squared = distance_squared(p, q);
                    sph_float w_r = smoothing_density(r_dist_squared);
                    p.density += q.mass * (w_r);
            }
        }
    std::cout<<"Particles density "<<p.density<<std::endl;

    }
}

void calculate_force_buckets(std::vector<particle> &particles, std::vector<std::vector<int>> &buckets, int nBuckets, sph_float dt)
{
    // Force = gravity_force + pressure_force + viscosity_force
    for (int i = 0; i < particles.size(); i++)
    {
        particle &p = particles[i];
        if (p.is_boundary)
            continue;
        int p_bucket_id = p.bucket_id;
        std::vector<int> neighbour_buckets = get_neighbours(p_bucket_id, nBuckets);
        for (int j = 0; j < neighbour_buckets.size(); j++)
        {
            int neighbour_bucket_id = neighbour_buckets[j];
            if (neighbour_bucket_id < 0 || neighbour_bucket_id >= nBuckets * nBuckets * nBuckets)
                continue;
            for (int k = 0; k < buckets[neighbour_bucket_id].size(); k++)
            {
                particle &q = particles[buckets[neighbour_bucket_id][k]];
                if (p.id < q.id || q.is_boundary) {
    pressure_viscosity_force(p, q, dt);
}

            }
        }
    }
}
