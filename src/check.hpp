

#include <vector>
#include <iostream>

using sph_float = float;
float PI = 3.14;
float r_e = 1.2;

struct particle
{
    int id;
    float x, y;
    float radius;
    float mass;
    float pressure;
    float density;
    float velocity_x, velocity_y;
    float force_x, force_y;
    int bucket_id = -1;

};

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
            p.id = i * cols + j;
            p.bucket_id = -1;
            particles.push_back(p);
        }
    }
    return particles;
}

std::vector<particle> particles = create_particles_grid(
    10, 10, // rows, cols
    0.05f,  // spacing between particles
    0.01f,  // radius of each particle
    1.0f    // mass of each particle
);

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

void apply_gravity(std::vector<particle> &particles, float gravity, float dt)
{
    for (auto &p : particles)
    {
        p.velocity_y += gravity * dt; // gravity affects only y-direction
    }
}

void update_velocities(std::vector<particle> &particles, float dt)
{
    for (auto &p : particles)
    {
        float acc_x = p.force_x / p.density;
        float acc_y = p.force_y / p.density;

        p.velocity_x += acc_x * dt;
        p.velocity_y += acc_y * dt;
    }
}

void update_positions(std::vector<particle> &particles, float dt)
{
    for (auto &p : particles)
    {
        p.x += p.velocity_x * dt;
        p.y += p.velocity_y * dt;
    }
}
void apply_boundary_conditions(std::vector<particle> &particles, float window_width, float window_height, float bounce = -0.8f)
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
    if (r < r_e)
    {
        sph_float val = 315 * std::pow((r_e * r_e - r), 3) / (64 * PI * std::pow(r_e, 9));
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

void calculate_pressure(std::vector<particle> &particles, sph_float k = 3, sph_float rho_o = 998.29, sph_float rest_pressure = 10.0)
{
    for (auto &p : particles)
    {
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
        if (r == 0.0f)
        {
            return 45.0f / (PI * std::pow(r_e, 6)); // tuned constant
        }
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
    // sph_float force_viscosity_x = q.mass*(q.vel_x - p.vel_x) * w_r_viscosity /(q.density);
    // sph_float force_viscosity_y = q.mass*(q.vel_y - p.vel_y) * w_r_viscosity /(q.density);
    sph_float fx = -force_pressure_x;
    sph_float fy = -force_pressure_y;

    p.force_x += fx;
    p.force_y += fy;

    q.force_x -= fx;
    q.force_y -= fy;
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
    int y = (idx / n) % n;
    int z = idx / (n * n);

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
                // sph_float dx = p.pos_x - q.pos_x;
                // sph_float dy = p.pos_y - q.pos_y;
                // sph_float dz = p.pos_z - q.pos_z;
                sph_float r_dist_squared = distance_squared(p, q);
                sph_float w_r = smoothing_density(r_dist_squared);
                p.density += q.mass * (w_r);
            }
        }
    }
}

void calculate_force_buckets(std::vector<particle> &particles, std::vector<std::vector<int>> &buckets, int nBuckets, sph_float dt)
{
    // Force = gravity_force + pressure_force + viscosity_force
    for(auto& p:particles){
        p.force_x =0;
        p.force_y=0;
    }
    for (int i = 0; i < particles.size(); i++)
    {
        particle &p = particles[i];
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
                if (p.id < q.id)
                    pressure_viscosity_force(p, q, dt);
            }
        }
    }
}

// void apply_repulsive_boundary_force(
//     std::vector<particle>& particles,
//     const std::vector<std::vector<int>>& neighbor_buckets,
//     int nBuckets
// ) {
//     int N = particles.size();

//     for (int a = 0; a < N; ++a) {
//         if (!particles[a].fix) {
//             particle& pa = particles[a];
//             sph_float dist = distance_squared(pa,pb);
//                     if (dist<r_e*r_e) { // boundary particle

//                         float dx = pa.x - pb.y;
//                         float dy = pa.y - pb.y;

//                         float r2 = dx*dx + dy*dy;
//                         if (r2 > r_e * r_e) continue;

//                         float r = std::sqrt(r2);
//                         if (r < 1e-5f) r = 1e-5f;

//                         // === Explicit spiky gradient ===
//                         float gradW_scalar = 45.0f * std::pow(r_e - r, 2) / (PI * std::pow(r_e, 6));
//                         float inv_r = 1.0f / r;

//                         float rho_a = pa.density;
//                         float rho_b = pb.density;
//                         float m_b = pb.mass;

//                         float term_a = (rho_a - rho0) / (rho_a * rho_a);
//                         float term_b = (rho_b - rho0) / (rho_b * rho_b);
//                         float coeff = -m_b * c0 * c0 * (term_a + term_b) * gradW_scalar * inv_r;

//                         sph_float force_x = coeff * dx;
//                         sph_float force_y = coeff * dy;

//                         pa.velocity_x += (-force_x/pa.density);
//                         pa.velocity_y += (-force_y/pa.density);
//                     }
//         }
//     }
// }

void apply_boundary_force(
    std::vector<particle>& particles,
    float window_width,
    float window_height)
{
    for (auto& p : particles) {
        // Check all 4 walls
        std::vector<std::pair<float, float>> wall_points;

        // Left wall
        if (p.x - p.radius < -window_width / 2.0f + r_e)
            wall_points.emplace_back(-window_width / 2.0f, p.y);

        // Right wall
        if (p.x + p.radius > window_width / 2.0f - r_e)
            wall_points.emplace_back(window_width / 2.0f, p.y);

        // Bottom wall
        if (p.y - p.radius < -window_height / 2.0f + r_e)
            wall_points.emplace_back(p.x, -window_height / 2.0f);

        // Top wall
        if (p.y + p.radius > window_height / 2.0f - r_e)
            wall_points.emplace_back(p.x, window_height / 2.0f);

        for (const auto& [bx, by] : wall_points) {
            float dx = p.x - bx;
            float dy = p.y - by;
            float r2 = dx * dx + dy * dy;

            if (r2 > r_e * r_e) continue;

            float r = std::sqrt(r2);
            if (r < 1e-5f) r = 1e-5f;

            // Spiky kernel gradient scalar part
            float gradW_scalar = 45.0f * std::pow(r_e - r, 2) / (PI * std::pow(r_e, 6));
            float inv_r = 1.0f / r;

            // Boundary repulsion force
            float rho = p.density;
            float m_b = p.mass; // or use fixed m_b = rho0 * V if needed

            float term = (rho - rho0) / (rho * rho);
            float coeff = -m_b * c0 * c0 * term * gradW_scalar * inv_r;

            float force_x = coeff * dx;
            float force_y = coeff * dy;

            p.velocity_x += (-force_x / p.density);
            p.velocity_y += (-force_y / p.density);
        }
    }
}