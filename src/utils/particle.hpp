#pragma once
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