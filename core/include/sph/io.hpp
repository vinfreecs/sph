#pragma once
// Legacy-VTK output (open in ParaView) and the VTK scene format written by tools/make_scene.py.

#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "particle.hpp"

namespace sph {

inline std::string frame_filename(int frame) { return "frame_" + std::to_string(frame) + ".vtk"; }

// Writes fluid particles only (ghost particles are static and would dominate the file).
inline bool write_vtk(const std::string& path, const std::vector<Particle>& ps) {
    std::ofstream out(path);
    if (!out) {
        std::cerr << "Failed to open " << path << "\n";
        return false;
    }
    size_t n = 0;
    for (const auto& p : ps) n += !p.fixed;

    out << "# vtk DataFile Version 4.0\nSPH particle data\nASCII\nDATASET UNSTRUCTURED_GRID\n";
    out << "POINTS " << n << " double\n" << std::fixed << std::setprecision(6);
    for (const auto& p : ps)
        if (!p.fixed) out << p.pos_x << " " << p.pos_y << " " << p.pos_z << "\n";
    out << "CELLS 0 0\nCELL_TYPES 0\nPOINT_DATA " << n << "\n";

    auto scalar = [&](const char* name, auto get) {
        out << "SCALARS " << name << " double\nLOOKUP_TABLE default\n";
        for (const auto& p : ps)
            if (!p.fixed) out << get(p) << "\n";
    };
    auto vector = [&](const char* name, auto get) {
        out << "VECTORS " << name << " double\n";
        for (const auto& p : ps)
            if (!p.fixed) get(p);
    };
    scalar("radius", [](const Particle& p) { return p.radius; });
    scalar("mass", [](const Particle& p) { return p.mass; });
    scalar("density", [](const Particle& p) { return p.density; });
    scalar("pressure", [](const Particle& p) { return p.pressure; });
    vector("velocity", [&](const Particle& p) { out << p.vel_x << " " << p.vel_y << " " << p.vel_z << "\n"; });
    vector("force", [&](const Particle& p) { out << p.force_x << " " << p.force_y << " " << p.force_z << "\n"; });
    return true;
}

// Reads POINTS + mass / velocity / radius / fixed arrays (the make_scene.py format).
inline bool read_vtk_scene(const std::string& path, std::vector<Particle>& ps) {
    std::ifstream in(path);
    if (!in) {
        std::cerr << "Failed to open " << path << "\n";
        return false;
    }
    std::string line;
    int n = 0;
    while (std::getline(in, line))
        if (line.rfind("POINTS", 0) == 0) {
            std::istringstream(line.substr(6)) >> n;
            break;
        }
    if (n <= 0) {
        std::cerr << path << ": no POINTS section\n";
        return false;
    }
    ps.assign(n, Particle{});
    for (int i = 0; i < n; ++i) {
        in >> ps[i].pos_x >> ps[i].pos_y >> ps[i].pos_z;
        ps[i].id = i;
    }
    auto seek = [&](const char* key, bool skip_lookup) {
        while (std::getline(in, line))
            if (line.find(key) != std::string::npos) {
                if (skip_lookup) std::getline(in, line);
                return true;
            }
        std::cerr << path << ": missing " << key << "\n";
        return false;
    };
    if (!seek("SCALARS mass", true)) return false;
    for (auto& p : ps) in >> p.mass;
    if (!seek("VECTORS velocity", false)) return false;
    for (auto& p : ps) in >> p.vel_x >> p.vel_y >> p.vel_z;
    if (!seek("SCALARS radius", true)) return false;
    for (auto& p : ps) in >> p.radius;
    if (!seek("SCALARS fixed", true)) return false;
    for (auto& p : ps) {
        int f = 0;
        in >> f;
        p.fixed = f != 0;
    }
    return static_cast<bool>(in);
}

}  // namespace sph
