__device__ sph_float distance_squared(particle p, particle q)
{
    sph_float dx = p.pos_x - q.pos_x;
    sph_float dy = p.pos_y - q.pos_y;
    sph_float dz = p.pos_z - q.pos_z;
    return dx * dx + dy * dy + dz * dz;
}

__device__ sph_float smoothing_pressure(sph_float r)
{
    if (r < r_e)
    {
        return 45.0f * powf(r_e - r, 3.0f) / (PI * powf(r_e, 6.0f) * r);
    }
    return 0.0f;
}

__device__ sph_float smoothing_viscosity(sph_float r)
{
    if (r < r_e)
    {
        return 45.0f * (r_e - r) / (PI * powf(r_e, 6.0f));
    }
    return 0.0f;
}

__device__ sph_float smoothing_density(sph_float r_squared)
{
    sph_float r = sqrtf(r_squared);
    if (r < r_e)
    {
        return 315.0f * powf((r_e * r_e - r * r), 3.0f) / (64.0f * PI * powf(r_e, 9.0f));
    }
    return 0.0f;
}

__global__ void reset_forces(particle *particles, int n)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i < n)
    {
        if(!particles[i].fixed){

            particles[i].force_x = 0;
            particles[i].force_y = 0;
            particles[i].force_z = 0;
        }
    }
}

__device__ void pressure_viscosity_force(particle &p, const particle &q)
{
    
    sph_float r_squared = distance_squared(p, q);
    sph_float r = sqrtf(r_squared);
    if (r > r_e)
        return;
    if(r < 1e-5f)
        r=1e-5f;
    sph_float w_r_pressure = smoothing_pressure(r);
    sph_float w_r_viscosity = smoothing_viscosity(r);

    sph_float dx = q.pos_x - p.pos_x;
    sph_float dy = q.pos_y - p.pos_y;
    sph_float dz = q.pos_z - p.pos_z;

    sph_float force_pressure_x = q.mass * (p.pressure + q.pressure) * w_r_pressure * dx / (2.0f * q.density);
    sph_float force_pressure_y = q.mass * (p.pressure + q.pressure) * w_r_pressure * dy / (2.0f * q.density);
    sph_float force_pressure_z = q.mass * (p.pressure + q.pressure) * w_r_pressure * dz / (2.0f * q.density);

    sph_float force_viscosity_x = q.mass * (q.vel_x - p.vel_x) * w_r_viscosity / q.density;
    sph_float force_viscosity_y = q.mass * (q.vel_y - p.vel_y) * w_r_viscosity / q.density;
    sph_float force_viscosity_z = q.mass * (q.vel_z - p.vel_z) * w_r_viscosity / q.density;

    p.force_x += -force_pressure_x + dynamic_viscosity * force_viscosity_x;
    p.force_y += -force_pressure_y + dynamic_viscosity * force_viscosity_y;
    p.force_z += -force_pressure_z + dynamic_viscosity * force_viscosity_z;
}

__device__ int get_hash(int cell_x, int cell_y, int cell_z)
{
    return (cell_x * 9283711) ^ (cell_y * 689287499) ^ (cell_z * 283923481);
}

__global__ void update_neighbors_kernel(
    particle *particles,
    int *neighbours,
    int *particle_neighbours,
    int num_particles)
{
    // the intialsization of neighbours and particle_neig should be done in a different kernal or use
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= num_particles)
        return;
    particle p = particles[i];
    int x = (int)((p.pos_x + boundary * 0.5f) / r_e);
    int y = (int)((p.pos_y + boundary * 0.5f) / r_e);
    int z = (int)((p.pos_z + boundary * 0.5f) / r_e);

    int hash = get_hash(x, y, z);
    int cell_id = ((hash % num_particles) + num_particles) % num_particles;

    // Atomically update hash table and get the previous value
    int previous = atomicExch(&neighbours[cell_id], i);
    particle_neighbours[i] = previous;
}

__global__ void compute_density_pressure_hash_kernel(
    particle *particles,
    const int *neighbours,
    const int *particle_neighbours,
    int num_particles)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= num_particles)
        return;

    particle &p = particles[i];
    p.density = rho_o;
    if(!p.fixed){
        
        p.density = 0.0f;
        p.pressure = 0.0f;
    
        int x = (int)((p.pos_x + boundary * 0.5f) / r_e);
        int y = (int)((p.pos_y + boundary * 0.5f) / r_e);
        int z = (int)((p.pos_z + boundary * 0.5f) / r_e);
    
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
                    int cell_id = ((hash % num_particles) + num_particles) % num_particles;
    
                    int j = neighbours[cell_id];
                    while (j != -1)
                    {
                        particle &q = particles[j];
                        sph_float r_sq = distance_squared(p, q);
                        sph_float w_r = smoothing_density(r_sq);
                        p.density += q.mass * w_r;
    
                        j = particle_neighbours[j];
                    }
                }
            }
        }
    }

    p.pressure = rest_pressure + k * (p.density - rho_o);
}

__global__ void calculate_force_hash_kernel(
    particle *particles,
    const int *neighbours,
    const int *particle_neighbours,
    int num_particles)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= num_particles)
        return;

    particle &p = particles[i];
    if(!p.fixed){
        p.force_x = 0;
        p.force_y = 0;
        p.force_z = 0;
    
        int x = (int)((p.pos_x + boundary * 0.5f) / r_e);
        int y = (int)((p.pos_y + boundary * 0.5f) / r_e);
        int z = (int)((p.pos_z + boundary * 0.5f) / r_e);
    
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
                    int cell_id = ((hash % num_particles) + num_particles) % num_particles;
    
                    int j = neighbours[cell_id];
                    while (j != -1)
                    {
                        if (j != i)
                        {
                            pressure_viscosity_force(p, particles[j]);
                        }
                        j = particle_neighbours[j];
                    }
                }
            }
        }

    }

    p.force_y -= p.density * gravity;
}



__global__ void update_velocity_kernel(particle* particles, int num_particles, sph_float dt)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= num_particles) return;

    particle& p = particles[i];
    if (p.density < 1e-5f) return;

    if(!p.fixed){
        p.vel_x += p.force_x * dt / p.density;
        p.vel_y += p.force_y * dt / p.density;
        p.vel_z += p.force_z * dt / p.density;

    }
}

__global__ void update_positions_kernel(particle* particles, int num_particles, sph_float dt)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= num_particles) return;

    particle& p = particles[i];

    if(!p.fixed){
    p.pos_x += p.vel_x * dt;
    p.pos_y += p.vel_y * dt;
    p.pos_z += p.vel_z * dt;}

}
