void apply_repulsive_boundary_force(
    std::vector<Particle>& particles,
    const std::vector<std::vector<int>>& neighbors
) {
    int N = particles.size();

    for (int a = 0; a < N; ++a) {
        if (particles[a].fix) continue; // Skip boundary particles

        // Loop over neighbors of a
        for (int b : neighbors[a]) {
            if (!particles[b].fix) continue; // Only act on boundary particles

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
