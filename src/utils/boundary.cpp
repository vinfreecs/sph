#include <vector>
#include <cmath>

//======================//
// Global Constants
//======================//
const double c0   = 20.0;          // Artificial speed of sound [m/s]
const double rho0 = 1000.0;        // Reference density [kg/m^3]
const double re   = 0.1;           // Kernel support radius

//======================//
// Basic 3D Vector Struct
//======================//
struct Vec3 {
    double x, y, z;
};

//======================//
// Particle Data Struct
//======================//
struct Particle {
    double x, y, z;      // Position
    double vx, vy, vz;   // Velocity
    double fx, fy, fz;   // Force accumulator
    double mass;         // Mass
    bool fix;            // true = boundary, false = fluid
    double density;      // Density ρ
    double pressure;     // Pressure P
};

//=========================================================//
// EOS Function: Linear Pressure-Density Relationship
//=========================================================//
// P = p0 + k*(rho - rho0) with k = c0^2
double compute_pressure(double rho) {
    const double p0 = 0.0;           // Optional background pressure
    const double k = c0 * c0;        // Stiffness coefficient

    return p0 + k * (rho - rho0);
}

//=========================================================//
// Kernel Gradient: ∇W_ab using Müller’s pressure kernel
//=========================================================//
// Computes ∇W(r_ab) where r_ab = a.pos - b.pos
Vec3 gradW_ab(const Particle& a, const Particle& b) {
    Vec3 r;
    r.x = a.x - b.x;
    r.y = a.y - b.y;
    r.z = a.z - b.z;

    double r_mag = std::sqrt(r.x * r.x + r.y * r.y + r.z * r.z);

    Vec3 grad = {0.0, 0.0, 0.0};

    if (r_mag > 0 && r_mag < re) {
        double coeff = 45.0 / (M_PI * std::pow(re, 6));
        double scale = coeff * std::pow(re - r_mag, 2) / r_mag;

        grad.x = scale * r.x;
        grad.y = scale * r.y;
        grad.z = scale * r.z;
    }

    return grad;
}

//=========================================================//
// Boundary Repulsion Force (Monaghan-style)
//=========================================================//
// Applies repulsive force to fluid particles from boundary particles
void apply_repulsive_boundary_force(std::vector<Particle>& particles) {
    int N = particles.size();

    for (int a = 0; a < N; ++a) {
        if (particles[a].fix) continue; // Skip boundary particles

        for (int b = 0; b < N; ++b) {
            if (!particles[b].fix) continue; // Only apply to boundary particles

            Vec3 gradW = gradW_ab(particles[a], particles[b]);

            double rho_a = particles[a].density;
            double rho_b = particles[b].density;
            double m_b   = particles[b].mass;

            double term_a = (rho_a - rho0) / (rho_a * rho_a);
            double term_b = (rho_b - rho0) / (rho_b * rho_b);

            double coeff = -m_b * c0 * c0 * (term_a + term_b);

            particles[a].fx += coeff * gradW.x;
            particles[a].fy += coeff * gradW.y;
            particles[a].fz += coeff * gradW.z;
        }
    }
}

//=========================================================//
// Pressure Update Wrapper
//=========================================================//
// Loops through particles and updates pressure using EOS
void update_pressure_from_density(std::vector<Particle>& particles) {
    for (auto& p : particles) {
        if (!p.fix) {
            p.pressure = compute_pressure(p.density);
        }
    }
}
