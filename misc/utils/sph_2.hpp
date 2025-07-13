// #pragma once
// #include "particle.hpp"
// #include "spatial_hashing.hpp"
// #include "dec.hpp"

// static const int neighbor_offsets[27][3] = {
//     // Bottom layer (z = -1)
//     {-1, -1, -1}, {-1, -1,  0}, {-1, -1,  1},
//     {-1,  0, -1}, {-1,  0,  0}, {-1,  0,  1},
//     {-1,  1, -1}, {-1,  1,  0}, {-1,  1,  1},
    
//     // Middle layer (z = 0)
//     { 0, -1, -1}, { 0, -1,  0}, { 0, -1,  1},
//     { 0,  0, -1}, { 0,  0,  0}, { 0,  0,  1},
//     { 0,  1, -1}, { 0,  1,  0}, { 0,  1,  1},
    
//     // Top layer (z = 1)
//     { 1, -1, -1}, { 1, -1,  0}, { 1, -1,  1},
//     { 1,  0, -1}, { 1,  0,  0}, { 1,  0,  1},
//     { 1,  1, -1}, { 1,  1,  0}, { 1,  1,  1}
// };

// sph_float distance_squared(particle p, particle q){
//     sph_float dx = p.pos_x -q.pos_x;
//     sph_float dy = p.pos_y -q.pos_y;
//     sph_float dz = p.pos_z -q.pos_z;
//     sph_float r_dist = dx*dx + dy*dy + dz*dz;
//     return r_dist;
// }

// sph_float smoothing_pressure(sph_float r){
//     if(r < r_e ){
//         sph_float val = 45 * std::pow(r_e - r,3)/(PI * std::pow(r_e,6)* r);
//         return val;
//     }
//     else return 0.0;
// }

// sph_float smoothing_viscosity(sph_float r){
//     if(r < r_e){
//         sph_float val = 45 * (r_e-r)/(PI * std::pow(r_e,6));
//         return val;
//     }
//     else return 0.0;
// }

// sph_float smoothing_density(sph_float r){
//     r = std::sqrt(r);
//     if(r<r_e){
//         sph_float val = 315*std::pow((r_e*r_e - r*r),3)/(64 * PI * std::pow(r_e,9));
//         return val;
//     }
//     else return 0.0;
// }

// void boundary_navie(particle &p, double boundary=3, double bounce_factor = 0.25)
// {
//     double hole_size = 0.5;
//     double half_boundary = boundary / 2.0;
//     double half_hole = hole_size / 2.0;

//     // X-axis
//     if (p.pos_x - radius < -half_boundary)
//     {
//         p.pos_x = -half_boundary + radius;
//         p.vel_x *= -bounce_factor;
//     }

//     else if (p.pos_x + radius > half_boundary)
//     {
//         p.pos_x = half_boundary - radius;
//         p.vel_x *= -bounce_factor;
//     }

//     // Y-axis
//     if (p.pos_y - radius < -half_boundary)
//     {
//         p.pos_y = -half_boundary + radius;
//         p.vel_y *= -bounce_factor;
//     }
//     else if (p.pos_y + radius > half_boundary)
//     {
//         p.pos_y = half_boundary - radius;
//         p.vel_y *= -bounce_factor;
//     }

//     // Z-axis
//     if (p.pos_z - radius < -half_boundary)
//     {
//         p.pos_z = -half_boundary + radius;
//         p.vel_z *= -bounce_factor;
//     }
//     else if (p.pos_z + radius > half_boundary)
//     {
//         p.pos_z = half_boundary - radius;
//         p.vel_z *= -bounce_factor;
//     }
// }

// void pressure_viscosity_force(particle& p, particle& q){
//     sph_float r_squared = distance_squared(p,q);
//     sph_float r = std::sqrt(r_squared);
//     if(r > r_e) return; // no interaction if distance is greater than smoothing length
//     if(r < 1e-5) r = 1e-5; // to avoid division by zero
//     // Calculate the force due to pressure and viscosity

//     sph_float w_r_pressure = smoothing_pressure(r);
//     sph_float w_r_viscosity = smoothing_viscosity(r);
//     sph_float force_pressure_x = q.mass*(p.pressure + q.pressure)* w_r_pressure * (q.pos_x -p.pos_x) / (2* q.density);//TODO
//     sph_float force_pressure_y = q.mass*(p.pressure + q.pressure)* w_r_pressure * (q.pos_y -p.pos_y) / (2* q.density);
//     sph_float force_pressure_z = q.mass*(p.pressure + q.pressure)* w_r_pressure * (q.pos_z -p.pos_z) / (2* q.density);
//     sph_float force_viscosity_x = q.mass*(q.vel_x - p.vel_x) * w_r_viscosity /(q.density);
//     sph_float force_viscosity_y = q.mass*(q.vel_y - p.vel_y) * w_r_viscosity /(q.density);
//     sph_float force_viscosity_z = q.mass*(q.vel_z - p.vel_z) * w_r_viscosity /(q.density);
//     sph_float fx = -force_pressure_x + dynamic_viscosity*force_viscosity_x;
//     sph_float fy= -force_pressure_y + dynamic_viscosity*force_viscosity_y;       
//     sph_float fz= -force_pressure_z + dynamic_viscosity*force_viscosity_z;   
//     p.force_x += fx;
//     p.force_y += fy;
//     p.force_z += fz;
//     q.force_x -= fx; // Newton's third law
//     q.force_y -= fy; // Newton's third law
//     q.force_z -= fz; // Newton's third law        
// }

// void update_velocity(std::vector<particle>& particles, sph_float dt){
//     for(auto &p : particles){
//         //std::cout <<"Force x :"<<p.force_x<<" y :"<<p.force_y<<" z :"<<p.force_z<<std::endl;
//         p.vel_x += p.force_x*dt/p.density;
//         p.vel_y += p.force_y*dt/p.density;
//         p.vel_z += ((p.force_z/p.density))*dt;

//     }
// }

// void update_positions(std::vector<particle>& particles, sph_float dt){
//     for(auto &p : particles){
//         p.pos_x += p.vel_x*dt;
//         p.pos_y += p.vel_y*dt;
//         p.pos_z += p.vel_z*dt;
//         boundary_navie(p);
//     }
// }

// void calculate_pressure(std::vector<particle>& particles, sph_float k =3, sph_float rho_o = 998.29, sph_float rest_pressure = 10.0){
//     for(auto & p: particles){
//         p.pressure = rest_pressure + k * (p.density - rho_o);
//         //std::cout<<"Particle "<<p.id<<" density :"<<p.density<<" pressure :"<<rest_pressure + k * (p.density - rho_o)<<std::endl;
//     }
// }


// std::vector<int> get_neighbours_hash(const particle& p, int table_size) {
//     int x = (int)((p.pos_x + (sph_float)boundary / 2.0f) / r_e);
//     int y = (int)((p.pos_y + (sph_float)boundary / 2.0f) / r_e);
//     int z = (int)((p.pos_z + (sph_float)boundary / 2.0f) / r_e);

//     std::vector<int> neighbor_hashes;

//     for (int dx = -1; dx <= 1; dx++) {
//         for (int dy = -1; dy <= 1; dy++) {
//             for (int dz = -1; dz <= 1; dz++) {
//                 int nx = x + dx;
//                 int ny = y + dy;
//                 int nz = z + dz;

//                 if (nx >= 0 && ny >= 0 && nz >= 0) {
//                     int hash = get_hash(nx, ny, nz);
//                     hash = (hash % table_size + table_size) % table_size;
//                     neighbor_hashes.push_back(hash);
//                 }
//             }
//         }
//     }

//     return neighbor_hashes;
// }

// std::vector<int> get_neighbours_hash_hardcoded(const particle& p, int table_size) {
//     const sph_float offset = boundary * 0.5f;
//     const sph_float inv_r_e = 1.0f / r_e;
    
//     const int base_x = (int)((p.pos_x + offset) * inv_r_e);
//     const int base_y = (int)((p.pos_y + offset) * inv_r_e);
//     const int base_z = (int)((p.pos_z + offset) * inv_r_e);

//     std::vector<int> neighbor_hashes;
//     neighbor_hashes.reserve(27);

//     for (int i = 0; i < 27; i++) {
//         const int nx = base_x + neighbor_offsets[i][0];
//         const int ny = base_y + neighbor_offsets[i][1];
//         const int nz = base_z + neighbor_offsets[i][2];

//         if (nx >= 0 && ny >= 0 && nz >= 0) {
//             int hash = get_hash(nx, ny, nz);
//             hash = ((hash % table_size) + table_size) % table_size;
//             neighbor_hashes.push_back(hash);
//         }
//     }

//     return neighbor_hashes;
// }

// void calculate_density_hash(std::vector<particle>& particles,
//                              const std::vector<Neighbour>& neighbours,
//                              const std::vector<int>& start_idx)
// {
//     int n = particles.size();

//     for (auto& p : particles) {
//         p.density = 0.0;
//         std::vector<int> neighbor_cells = get_neighbours_hash_hardcoded(p, n);

//         for (int cell_idx : neighbor_cells) {
//             int start = start_idx[cell_idx];
//             if (start == -1) continue;

//             for (int k = start; k < neighbours.size() && neighbours[k].cell_id == cell_idx; ++k) {
//                 const particle& q = particles[neighbours[k].particle_id];
//                 sph_float r2 = distance_squared(p, q);
//                 sph_float w = smoothing_density(r2);
//                 p.density += q.mass * w;
//             }
//         }
//     }
// }


// void calculate_force_hash(std::vector<particle>& particles,
//                           const std::vector<Neighbour>& neighbours,
//                           const std::vector<int>& start_idx)
// {
//     int n = particles.size();

//     for (auto& p : particles) {
//         p.force_x = p.force_y = p.force_z = 0.0;
//     }

//     for (auto& p : particles) {
//         std::vector<int> neighbor_cells = get_neighbours_hash_hardcoded(p, n);

//         for (int cell_idx : neighbor_cells) {
//             int start = start_idx[cell_idx];
//             if (start == -1) continue;

//             for (int k = start; k < neighbours.size() && neighbours[k].cell_id == cell_idx; ++k) {
//                 particle& q = particles[neighbours[k].particle_id];
//                 if (p.id < q.id) pressure_viscosity_force(p, q);
//             }
//         }

//         // Apply gravity
//         p.force_z -= p.mass * gravity;
//     }
// }





