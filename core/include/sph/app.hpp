#pragma once
// Command line + scene setup shared by the CPU and CUDA executables.

#include <sys/stat.h>

#include <cerrno>
#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>

#include "io.hpp"
#include "params.hpp"
#include "scenes.hpp"

namespace sph {

struct Scene {
    int dim = 3;
    Params params;
    std::vector<Particle> particles;
    sph_float dt = 1e-4f;
    sph_float t_end = 5.0f;
};

struct RunConfig {
    std::string scene = "box";
    int fluid = 10000;
    long steps = -1;           // overrides t_end when >= 0
    sph_float dt = 0;          // overrides the scene's dt when > 0
    int every = 500;           // VTK output interval in steps, 0 disables output
    std::string out = "output";
};

inline void print_usage(const char* prog) {
    std::cout << "usage: " << prog << " [options]\n"
              << "  --scene NAME|FILE.vtk  box (3D, default), cup (2D), or a make_scene.py VTK file\n"
              << "  --fluid N              fluid particles for generated scenes (default 10000)\n"
              << "  --steps N              number of time steps (default: scene end time / dt)\n"
              << "  --dt DT                time step override\n"
              << "  --every N              write a VTK frame every N steps, 0 = never (default 500)\n"
              << "  --out DIR              output directory (default ./output)\n";
}

inline RunConfig parse_args(int argc, char** argv) {
    RunConfig c;
    for (int i = 1; i < argc; ++i) {
        std::string a = argv[i];
        auto val = [&]() -> std::string {
            if (i + 1 >= argc) {
                std::cerr << "missing value for " << a << "\n";
                std::exit(2);
            }
            return argv[++i];
        };
        if (a == "--scene") c.scene = val();
        else if (a == "--fluid") c.fluid = std::stoi(val());
        else if (a == "--steps") c.steps = std::stol(val());
        else if (a == "--dt") c.dt = std::stof(val());
        else if (a == "--every") c.every = std::stoi(val());
        else if (a == "--out") c.out = val();
        else if (a == "-h" || a == "--help") {
            print_usage(argv[0]);
            std::exit(0);
        } else {
            std::cerr << "unknown option " << a << "\n";
            print_usage(argv[0]);
            std::exit(2);
        }
    }
    return c;
}

// mkdir -p (POSIX). Avoids std::filesystem, which needs -lstdc++fs on older GCC.
inline bool make_dirs(const std::string& path) {
    for (size_t pos = 0; pos != std::string::npos;) {
        pos = path.find('/', pos + 1);
        std::string part = path.substr(0, pos);
        if (!part.empty() && mkdir(part.c_str(), 0755) != 0 && errno != EEXIST) {
            std::cerr << "cannot create directory " << part << "\n";
            return false;
        }
    }
    return true;
}

inline bool ends_with(const std::string& s, const std::string& suffix) {
    return s.size() >= suffix.size() && s.compare(s.size() - suffix.size(), suffix.size(), suffix) == 0;
}

// box  : the HESP submission setup — fluid cube (edge 0.5) inside a closed ghost box (edge 1.5).
// cup  : 2D fluid square falling onto a cup (ring with a 30 degree opening) inside a ghost box.
// *.vtk: scene from tools/make_scene.py (free fall / flow over sphere).
inline bool make_scene(const std::string& name, int fluid, Scene& s) {
    if (name == "box") {
        s.dim = 3;
        s.params = preset_box_3d();
        s.particles = fluid_cube_3d(fluid, 0.5, s.params.particle_radius);
        append(s.particles, boundary_box_3d(1.5f, s.params.boundary_spacing));
        s.dt = 1e-4f;
        s.t_end = 5.0f;
    } else if (name == "cup") {
        s.dim = 2;
        s.params = preset_cup_2d();
        s.particles = fluid_square_2d(fluid, 1.0, s.params.particle_radius);
        append(s.particles, boundary_box_2d(1.5f, s.params.boundary_spacing));
        append(s.particles, boundary_ring_2d(0.25f, 0.04f, 0.2f, kPi / 6, kPi / 2));
        s.dt = 6.25e-4f;
        s.t_end = 5.0f;
    } else if (ends_with(name, ".vtk")) {
        s.dim = 3;
        s.params = preset_vtk_3d();
        if (!read_vtk_scene(name, s.particles)) return false;
        s.dt = 1e-4f;
        s.t_end = 5.0f;
    } else {
        std::cerr << "unknown scene '" << name << "'\n";
        return false;
    }
    assign_ids(s.particles);
    return true;
}

inline void print_scene_summary(const Scene& s) {
    size_t fluid = 0;
    for (const auto& p : s.particles) fluid += !p.fixed;
    std::cout << s.dim << "D scene: " << fluid << " fluid + " << s.particles.size() - fluid
              << " boundary = " << s.particles.size() << " particles, dt = " << s.dt << "\n";
}

}  // namespace sph
