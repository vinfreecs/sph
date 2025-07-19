using sph_float = float; // just to quickly change from float to double

extern __constant__ sph_float r_e = 0.04; // assuming 6 times of the radius this is the smoothing radius
 sph_float r_e_b = 0.02; // the radius of boundary particles
extern __constant__ sph_float PI = 3.14;
extern __constant__ sph_float dynamic_viscosity = 10;
extern __constant__ sph_float gravity = 9.8;
extern __constant__ sph_float radius = 0.03;
extern __constant__ sph_float k = 3; // pressure multiplier
extern __constant__ sph_float rho_o = 998.29; 
extern __constant__ sph_float rest_pressure = 0.0;
extern __constant__ int boundary = 15;

sph_float time_step = 0.0001f;
sph_float total_time = 5.0f;

sph_float output_freq = 500.0f;

sph_float radius_device = 0.03; 


