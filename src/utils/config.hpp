using sph_float = float;

extern __constant__ sph_float r_e = 0.04; // assuming 6 times of the radius
extern __constant__ sph_float r_e_b = 0.02;
extern __constant__ sph_float PI = 3.14;
extern __constant__ sph_float dynamic_viscosity = 10;
extern __constant__ sph_float gravity = 9.8;
extern __constant__ sph_float radius = 0.03;
extern __constant__ sph_float k = 3;
extern __constant__ sph_float rho_o = 998.29;
extern __constant__ sph_float rest_pressure = 0.0;
extern __constant__ int boundary = 2;
