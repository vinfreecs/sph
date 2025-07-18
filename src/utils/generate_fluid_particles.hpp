#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <string>
#include <filesystem>


std::vector<particle> create_boundary_particles_3d_aks(sph_float boundary_size, sph_float r_e_b, float overlap_factor = 0.2f) {
    std::vector<particle> boundary_particles;

    sph_float spacing = r_e_b * overlap_factor;
    int num_per_axis = static_cast<int>(boundary_size / spacing);
    sph_float half_size = boundary_size / 2.0f;

    int id = 0;
    for (int i = 0; i <= num_per_axis; ++i) {
        for (int j = 0; j <= num_per_axis; ++j) {
            for (int k = 0; k <= num_per_axis; ++k) {
                if (k == num_per_axis) continue; // Skip top wall

                bool is_border = (
                    i == 0 || i == num_per_axis ||
                    j == 0 || j == num_per_axis ||
                    k == 0
                );

                if (is_border) {
                    sph_float x = -half_size + i * spacing;
                    sph_float y = -half_size + j * spacing;
                    sph_float z = -half_size + k * spacing;

                    particle p;
                    p.id = id++;
                    p.pos_x = x;
                    p.pos_y = y;
                    p.pos_z = z;
                    p.vel_x = p.vel_y = p.vel_z = 0.0f;
                    p.mass = 1.0f;
                    p.density = 7870.0f;
                    p.pressure = 0.0f;
                    p.force_x = p.force_y = p.force_z = 0.0f;
                    p.fixed = true;
                    p.radius = r_e_b / 2.0f;
                    boundary_particles.push_back(p);
                }
            }
        }
    }

    return boundary_particles;
}

std::vector<particle> create_fixed_fluid_block_aks(int per_axis, double spacing, double solid_cube_size) {
    std::vector<particle> particles_;

    double fluid_cube_size = per_axis * spacing;
    double half_fluid = fluid_cube_size / 2.0;

    // Offset fluid upward so half is outside solid
    double z_offset = solid_cube_size / 4.0;

    int id = 0;
    for (int i = 0; i < per_axis; ++i) {
        for (int j = 0; j < per_axis; ++j) {
            for (int k = 0; k < per_axis; ++k) {
                particle p;
                p.id = id++;
                p.pos_x = i * spacing + spacing / 2 - half_fluid;
                p.pos_y = j * spacing + spacing / 2 - half_fluid;
                p.pos_z = k * spacing + spacing / 2 - half_fluid + z_offset;
                p.vel_x = p.vel_y = p.vel_z = 0.0f;
                p.mass = 1.0f;
                p.density = 1000.0f;
                p.pressure = 0.0f;
                p.force_x = p.force_y = p.force_z = 0.0f;
                p.fixed = false;
                p.radius = spacing / 2;
                particles_.push_back(p);
            }
        }
    }

    std::cout << "✅ Created " << particles_.size() << " fluid particles (≈50k)\n";
    return particles_;
}

// int main() {
//     // Reference spacing
//     int ref_particles = 10000;
//     double ref_cube_size = 1.0;
//     int ref_axis = static_cast<int>(std::ceil(std::cbrt(ref_particles)));
//     double spacing = ref_cube_size / ref_axis;

//     // === SETTINGS ===
//     double solid_cube_size = 3.0;     // Bigger domain
//     int fluid_per_axis = 41;          // ~50,000 fluid particles

//     // === GENERATE SYSTEM ===
//     auto boundary = create_boundary_particles_3d(solid_cube_size, spacing);
//     auto fluid = create_fixed_fluid_block(fluid_per_axis, spacing, solid_cube_size);

//     std::cout << "Solid (boundary) particles: " << boundary.size() << "\n";
//     std::cout << "Fluid particles: " << fluid.size() << "\n";
//     std::cout << "Total particles: " << boundary.size() + fluid.size() << "\n";

//     std::vector<particle> all = fluid;
//     all.insert(all.end(), boundary.begin(), boundary.end());

//     write_particles_to_vtk("./output/particles_fixed_fluid_test.vtk", all);
//     return 0;
// }

