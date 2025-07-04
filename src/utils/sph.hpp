using sph_float = float;

sph_float r_e = 0.006; // assuming 6 times of the radius
sph_float PI = 3.14;
sph_float dynamic_viscosity = 0.001;
sph_float gravity = 9.8;
sph_float radius = 0.001;

sph_float distance_squared(particle p, particle q){
    sph_float dx = p.pos_x -q.pos_x;
    sph_float dy = p.pos_y -q.pos_y;
    sph_float dz = p.pos_z -q.pos_z;
    sph_float r_dist = dx*dx + dy+dy+dz*dz;
    return r_dist;
}

sph_float smoothing_pressure(sph_float r){
    if(r < r_e){
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
    if(r<r_e){
        sph_float val = 315*std::pow(r_e*r_e - r,3)/(64 * PI * std::pow(r_e,9));
        return val;
    }
    else return 0.0;
}

void boundary_navie(particle &p, double boundary=3, double bounce_factor = 0.5)
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


void calculate_density(std::vector<particle>& particles){
    for(auto & p : particles){
        p.density = 0.0;
        for(auto& q : particles){
            if(p.id != q.id){
                sph_float dx = p.pos_x -q.pos_x;
                sph_float dy = p.pos_y -q.pos_y;
                sph_float dz = p.pos_z -q.pos_z;
                sph_float r_dist_squared = distance_squared(p,q);
                sph_float w_r = smoothing_density(r_dist_squared);
                p.density+= q.mass * (w_r);
            }
        }
    }
}

void calculate_pressure(std::vector<particle>& particles, sph_float k =1000.0, sph_float rho_o = 1000.0, sph_float rest_pressure = 0.0){
    for(auto & p: particles){
        p.pressure = rest_pressure + k * (p.density - rho_o);
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
                sph_float r_squared = distance_squared(p,q);
                sph_float r = std::sqrt(r_squared);
                sph_float w_r_pressure = smoothing_pressure(r);
                sph_float w_r_viscosity = smoothing_viscosity(r);
                sph_float force_pressure_x = q.mass*(p.pressure + q.pressure)* w_r_pressure * (p.pos_x -q.pos_x) / (2* q.density);
                sph_float force_pressure_y = q.mass*(p.pressure + q.pressure)* w_r_pressure * (p.pos_y -q.pos_y) / (2* q.density);
                sph_float force_pressure_z = q.mass*(p.pressure + q.pressure)* w_r_pressure * (p.pos_z -q.pos_z) / (2* q.density);
                sph_float force_viscosity_x = q.mass*(q.vel_x - p.vel_x) * w_r_viscosity /(q.density);
                sph_float force_viscosity_y = q.mass*(q.vel_y - p.vel_y) * w_r_viscosity /(q.density);
                sph_float force_viscosity_z = q.mass*(q.vel_z - p.vel_z) * w_r_viscosity /(q.density);
                p.force_x += -force_pressure_x + dynamic_viscosity*force_viscosity_x;
                p.force_y += -force_pressure_y + dynamic_viscosity*force_viscosity_y;               
                p.force_z += -force_pressure_z + dynamic_viscosity*force_viscosity_z;               
            }
            
        }
        p.force_z -= p.mass * gravity; // gravity only in the -ve z direction
    }
}

void update_velocity(std::vector<particle>& particles, sph_float dt){
    for(auto &p : particles){
        //std::cout <<"Force x :"<<p.force_x<<" y :"<<p.force_y<<" z :"<<p.force_z<<std::endl;
        p.vel_x += p.force_x*dt/p.mass;
        p.vel_y += p.force_y*dt/p.mass;
        p.vel_z += p.force_z*dt/p.mass;

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

