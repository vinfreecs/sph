#include "dec.hpp"

void update_buckets(std::vector<particle>& particles, std::vector<std::vector<int>>& buckets, int boundary,int nBuckets){
    for(int i=0; i<particles.size();i++){
        sph_float x = particles[i].pos_x + (sph_float)boundary/2;
        sph_float y = particles[i].pos_y + (sph_float)boundary/2;
        sph_float z = particles[i].pos_z + (sph_float)boundary/2;
        if(x < 0 || x > boundary || y < 0 || y > boundary || z < 0 || z > boundary){
            std::cerr<<"Particle "<<i<<" is out of bounds. Please check the input file."<<std::endl;
            std::cerr<<"x: "<<x<<" y: "<<y<<" z: "<<z<<" bucket_idx: "<<std::endl;
            std::cerr<<"nBuckets: "<<nBuckets<<" boundary: "<<boundary<<std::endl;
            continue;
        }
        int bucket_idx = (int)(x/r_e) + int(y/r_e)*nBuckets + int(z/r_e)*nBuckets*nBuckets;
        
        if(bucket_idx > nBuckets*nBuckets*nBuckets-1 || bucket_idx < 0){
            std::cerr<<"Particle "<<i<<" is out of bounds. Please check the input file."<<std::endl;
            std::cerr<<"x: "<<x<<" y: "<<y<<" z: "<<z<<" bucket_idx: "<<bucket_idx<<std::endl;
            std::cerr<<"nBuckets: "<<nBuckets<<" boundary: "<<boundary<<std::endl;
            continue;
            
        }
        if(particles[i].bucket_id == bucket_idx){
            continue; // already in the same bucket
        }
        else{
            if(particles[i].bucket_id !=-1){
                // Removing from old bucket
                int old_bucket_idx = particles[i].bucket_id;
                for(int j = 0; j<buckets[old_bucket_idx].size();j++){
                    if(buckets[old_bucket_idx][j] == particles[i].id){
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

std::vector<int> get_neighbours(int idx, int n) {
    int x = idx % n;
    int y = (idx / n) % n;
    int z = idx / (n * n);

    std::vector<int> neighbors;
    for (int dx = -1; dx <= 1; dx++) {
        for (int dy = -1; dy <= 1; dy++) {
            for (int dz = -1; dz <= 1; dz++) {
                int nx = x + dx;
                int ny = y + dy;
                int nz = z + dz;
                if (nx >= 0 && ny >= 0 && nz >= 0 && nx < n && ny < n && nz < n) {
                    int neighbor_idx = nx + ny * n + nz * n * n;
                    neighbors.push_back(neighbor_idx);
                }
            }
        }
    }
    return neighbors;
}


sph_float distance_squared(particle p, particle q){
    sph_float dx = p.pos_x -q.pos_x;
    sph_float dy = p.pos_y -q.pos_y;
    sph_float dz = p.pos_z -q.pos_z;
    sph_float r_dist = dx*dx + dy*dy+dz*dz;
    return r_dist;
}

sph_float smoothing_pressure(sph_float r){
    if(r < 1e-5) r = 1e-5; // to avoid division by zero
    if(r < r_e ){
        sph_float val = 45 * std::pow(r_e - r,3)/(PI * std::pow(r_e,6)* r);
        return val;
    }
    else return 0.0;
}

sph_float smoothing_viscosity(sph_float r){
    if(r_e > r){
        sph_float val = 45 * (r_e-r)/(PI * std::pow(r_e,6));
        return val;
    }
    else return 0.0;
}

sph_float smoothing_density(sph_float r){
    r = std::sqrt(r);
    if(r<r_e){
        sph_float val = 315*std::pow((r_e*r_e - r*r),3)/(64 * PI * std::pow(r_e,9));
        return val;
    }
    else return 0.0;
}

void boundary_navie(particle &p, double boundary=3, double bounce_factor = 0.25)
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

void pressure_viscosity_force(particle& p, particle& q){
    sph_float r_squared = distance_squared(p,q);
    sph_float r = std::sqrt(r_squared);
    if(r > r_e) return; // no interaction if distance is greater than smoothing length
    if(r < 1e-5) r = 1e-5; // to avoid division by zero
    // Calculate the force due to pressure and viscosity

    sph_float w_r_pressure = smoothing_pressure(r);
    sph_float w_r_viscosity = smoothing_viscosity(r);
    sph_float force_pressure_x = q.mass*(p.pressure + q.pressure)* w_r_pressure * (q.pos_x -p.pos_x) / (2* q.density);//TODO
    sph_float force_pressure_y = q.mass*(p.pressure + q.pressure)* w_r_pressure * (q.pos_y -p.pos_y) / (2* q.density);
    sph_float force_pressure_z = q.mass*(p.pressure + q.pressure)* w_r_pressure * (q.pos_z -p.pos_z) / (2* q.density);
    sph_float force_viscosity_x = q.mass*(q.vel_x - p.vel_x) * w_r_viscosity /(q.density);
    sph_float force_viscosity_y = q.mass*(q.vel_y - p.vel_y) * w_r_viscosity /(q.density);
    sph_float force_viscosity_z = q.mass*(q.vel_z - p.vel_z) * w_r_viscosity /(q.density);
    sph_float fx = -force_pressure_x + dynamic_viscosity*force_viscosity_x;
    sph_float fy= -force_pressure_y + dynamic_viscosity*force_viscosity_y;       
    sph_float fz= -force_pressure_z + dynamic_viscosity*force_viscosity_z;   
    p.force_x += fx;
    p.force_y += fy;
    p.force_z += fz;
    q.force_x -= fx; // Newton's third law
    q.force_y -= fy; // Newton's third law
    q.force_z -= fz; // Newton's third law        
}


void calculate_density(std::vector<particle>& particles){
    for(auto & p : particles){
        p.density = 0.0;
        for(auto& q : particles){
            sph_float dx = p.pos_x - q.pos_x;
            sph_float dy = p.pos_y - q.pos_y;
            sph_float dz = p.pos_z - q.pos_z;
            sph_float r_dist_squared = distance_squared(p, q);
            sph_float w_r = smoothing_density(r_dist_squared);
            p.density += q.mass * (w_r);
        }
    }
}

void calculate_density_buckets(std::vector<particle>& particles, std::vector<std::vector<int>>& buckets,int nBuckets){
    for(auto& p :particles){
        p.density = 0.0;
        int p_bucket_id = p.bucket_id;
        std::vector<int> neighbour_buckets = get_neighbours(p_bucket_id,nBuckets);
        for(int j = 0; j<neighbour_buckets.size(); j++){
            int neighbour_bucket_id = neighbour_buckets[j];
            if(neighbour_bucket_id < 0 || neighbour_bucket_id >= nBuckets*nBuckets*nBuckets) continue;
            for(int k =0 ; k<buckets[neighbour_bucket_id].size();k++){
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

void calculate_pressure(std::vector<particle>& particles, sph_float k =3, sph_float rho_o = 998.29, sph_float rest_pressure = 10.0){
    for(auto & p: particles){
        p.pressure = rest_pressure + k * (p.density - rho_o);
        //std::cout<<"Particle "<<p.id<<" density :"<<p.density<<" pressure :"<<rest_pressure + k * (p.density - rho_o)<<std::endl;
    }
}

void calculate_force(std::vector<particle>& particles){
    // Force = gravity_force + pressure_force + viscosity_force
    for(auto & p : particles){
        p.force_x=0;
        p.force_y=0;
        p.force_z=0;
    }
    for(int i =0 ; i <particles.size(); i++){
        particle &p = particles[i];
        for(int j = 0; j<particles.size(); j++){
            particle &q = particles[j];
            if(p.id != q.id){
                pressure_viscosity_force(p,q);            
            }
            
        }
        p.force_z -= p.density * gravity; // gravity only in the -ve z direction
    }
}

// std::vector<int> get_neighbours(int p_bucket_id, int nBuckets){
//     std::vector<int> neighbour_buckets(27);
//     neighbour_buckets[0] = p_bucket_id;
//     neighbour_buckets[1] = p_bucket_id - 1;
//     neighbour_buckets[2] = p_bucket_id + 1;
//     neighbour_buckets[3] = p_bucket_id - nBuckets;
//     neighbour_buckets[4] = p_bucket_id + nBuckets;
//     neighbour_buckets[5] = p_bucket_id - nBuckets - 1;
//     neighbour_buckets[6] = p_bucket_id - nBuckets + 1;
//     neighbour_buckets[7] = p_bucket_id + nBuckets - 1;
//     neighbour_buckets[8] = p_bucket_id + nBuckets + 1;
//     neighbour_buckets[9] = p_bucket_id - nBuckets*nBuckets;
//     neighbour_buckets[10] = p_bucket_id + nBuckets*nBuckets;
//     neighbour_buckets[11] = p_bucket_id - nBuckets*nBuckets - 1;
//     neighbour_buckets[12] = p_bucket_id - nBuckets*nBuckets + 1;
//     neighbour_buckets[13] = p_bucket_id + nBuckets*nBuckets - 1;
//     neighbour_buckets[14] = p_bucket_id + nBuckets*nBuckets + 1;
//     neighbour_buckets[15] = p_bucket_id - nBuckets*nBuckets - nBuckets;
//     neighbour_buckets[16] = p_bucket_id - nBuckets*nBuckets + nBuckets;
//     neighbour_buckets[17] = p_bucket_id + nBuckets*nBuckets - nBuckets;
//     neighbour_buckets[18] = p_bucket_id + nBuckets*nBuckets + nBuckets;
//     neighbour_buckets[19] = p_bucket_id - nBuckets*nBuckets - nBuckets - 1;
//     neighbour_buckets[20] = p_bucket_id - nBuckets*nBuckets - nBuckets + 1;
//     neighbour_buckets[21] = p_bucket_id - nBuckets*nBuckets + nBuckets - 1;
//     neighbour_buckets[22] = p_bucket_id - nBuckets*nBuckets + nBuckets + 1;
//     neighbour_buckets[23] = p_bucket_id + nBuckets*nBuckets - nBuckets - 1;
//     neighbour_buckets[24] = p_bucket_id + nBuckets*nBuckets - nBuckets + 1;
//     neighbour_buckets[25] = p_bucket_id + nBuckets*nBuckets + nBuckets - 1;
//     neighbour_buckets[26] = p_bucket_id + nBuckets*nBuckets + nBuckets + 1;
//     return neighbour_buckets;
// }


void calculate_force_buckets(std::vector<particle>& particles, std::vector<std::vector<int>>& buckets, int nBuckets){
    // Force = gravity_force + pressure_force + viscosity_force
    for(auto & p : particles){
        p.force_x=0;
        p.force_y=0;
        p.force_z=0;
    }
    for(int i =0 ; i <particles.size(); i++){
        particle &p = particles[i];
        int p_bucket_id = p.bucket_id;
        std::vector<int> neighbour_buckets = get_neighbours(p_bucket_id, nBuckets);
        for(int j = 0; j<neighbour_buckets.size(); j++){
            int neighbour_bucket_id = neighbour_buckets[j];
            if(neighbour_bucket_id < 0 || neighbour_bucket_id >= nBuckets*nBuckets*nBuckets) continue;
            for(int k =0 ; k<buckets[neighbour_bucket_id].size();k++){
                particle &q = particles[buckets[neighbour_bucket_id][k]];
                if(p.id < q.id) pressure_viscosity_force(p,q);
            }
        }
        p.force_z -= p.density * gravity; // gravity only in the -ve z direction 
        //TODO
    }
}

void update_velocity(std::vector<particle>& particles, sph_float dt){
    for(auto &p : particles){
        //std::cout <<"Force x :"<<p.force_x<<" y :"<<p.force_y<<" z :"<<p.force_z<<std::endl;
        p.vel_x += p.force_x*dt/p.density;
        p.vel_y += p.force_y*dt/p.density;
        p.vel_z += ((p.force_z/p.density))*dt;

    }
}

void update_positions(std::vector<particle>& particles, sph_float dt){
    for(auto &p : particles){
        p.pos_x += p.vel_x*dt;
        p.pos_y += p.vel_y*dt;
        p.pos_z += p.vel_z*dt;
        boundary_navie(p);
    }
}

void apply_repulsive_boundary_force(
    std::vector<particle>& particles,
    const std::vector<std::vector<int>>& neighbor_buckets,
    int nBuckets
) {
    int N = particles.size();

    for (int a = 0; a < N; ++a) {
        if (!particles[a].fix) {
            particle& pa = particles[a];
            int bucket_id = pa.bucket_id;
            std::vector<int> neighbours = get_neighbours(bucket_id, nBuckets);

            for (int b_bucket : neighbours) {
                if (b_bucket < 0 || b_bucket >= nBuckets * nBuckets * nBuckets) continue;

                for (int b_index : neighbor_buckets[b_bucket]) {
                    particle& pb = particles[b_index];
                    sph_float dist = distance_squared(pa,pb);
                    if (dist<r_e*r_e) { // boundary particle

                        float dx = pa.pos_x - pb.pos_x;
                        float dy = pa.pos_y - pb.pos_y;
                        float dz = pa.pos_z - pb.pos_z;

                        float r2 = dx*dx + dy*dy + dz*dz;
                        if (r2 > r_e * r_e) continue;

                        float r = std::sqrt(r2);
                        if (r < 1e-5f) r = 1e-5f;

                        // === Explicit spiky gradient ===
                        float gradW_scalar = 45.0f * std::pow(r_e - r, 2) / (PI * std::pow(r_e, 6));
                        float inv_r = 1.0f / r;

                        float rho_a = pa.density;
                        float rho_b = pb.density;
                        float m_b = pb.mass;

                        float term_a = (rho_a - rho0) / (rho_a * rho_a);
                        float term_b = (rho_b - rho0) / (rho_b * rho_b);
                        float coeff = -m_b * c0 * c0 * (term_a + term_b) * gradW_scalar * inv_r;

                        pa.force_x += coeff * dx;
                        pa.force_y += coeff * dy;
                        pa.force_z += coeff * dz;
                    }
                }
            }
        }
    }
}
