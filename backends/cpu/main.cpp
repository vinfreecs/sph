// Serial CPU reference: same physics as the CUDA backend, one particle at a time.

#include <chrono>
#include <iostream>

#include "sph/app.hpp"
#include "sph/cpu_solver.hpp"

template <int D>
static void run(sph::Scene& s, const sph::RunConfig& cfg) {
    sph::CpuSolver<D> solver;
    const long steps = cfg.steps >= 0 ? cfg.steps : static_cast<long>(s.t_end / s.dt);
    double compute = 0;
    if (cfg.every > 0) sph::write_vtk(cfg.out + "/" + sph::frame_filename(0), s.particles);
    for (long t = 0; t < steps; ++t) {
        auto t0 = std::chrono::steady_clock::now();
        solver.step(s.particles, s.params, s.dt);
        compute += std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
        if (cfg.every > 0 && (t + 1) % cfg.every == 0) {
            sph::write_vtk(cfg.out + "/" + sph::frame_filename(int(t + 1)), s.particles);
            std::cout << "step " << t + 1 << "/" << steps << "  compute " << compute << " s\n";
        }
    }
    std::cout << "Simulation time (excluding I/O): " << compute << " s for " << steps << " steps\n";
}

int main(int argc, char** argv) {
    sph::RunConfig cfg = sph::parse_args(argc, argv);
    sph::Scene scene;
    if (!sph::make_scene(cfg.scene, cfg.fluid, scene)) return 1;
    if (cfg.dt > 0) scene.dt = cfg.dt;
    if (cfg.every > 0 && !sph::make_dirs(cfg.out)) return 1;
    sph::print_scene_summary(scene);
    if (scene.dim == 2) run<2>(scene, cfg);
    else run<3>(scene, cfg);
}
