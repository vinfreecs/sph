// CUDA backend: one thread per particle, cell-linked list built with atomicExch.
// Physics lives in core/include/sph/physics.hpp and is shared with the CPU backend.

#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <vector>

#include <cuda_runtime.h>

#include "sph/app.hpp"
#include "sph/physics.hpp"

#define CUDA_CHECK(call)                                                                       \
    do {                                                                                       \
        cudaError_t err_ = (call);                                                             \
        if (err_ != cudaSuccess) {                                                             \
            std::fprintf(stderr, "CUDA error %s:%d: %s\n", __FILE__, __LINE__,                 \
                         cudaGetErrorString(err_));                                            \
            std::exit(1);                                                                      \
        }                                                                                      \
    } while (0)

using sph::Params;
using sph::Particle;
using sph::sph_float;

template <int D>
__global__ void build_cells_kernel(const Particle* ps, int* head, int* next, int n, sph_float h) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= n) return;
    int slot = sph::home_slot<D>(ps[i], h, n);
    next[i] = atomicExch(&head[slot], i);
}

template <int D>
__global__ void density_pressure_kernel(Particle* ps, const int* head, const int* next, int n, Params P) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i < n) sph::compute_density_pressure<D>(i, ps, head, next, n, P);
}

template <int D>
__global__ void force_kernel(Particle* ps, const int* head, const int* next, int n, Params P) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i < n) sph::compute_forces<D>(i, ps, head, next, n, P);
}

template <int D>
__global__ void integrate_kernel(Particle* ps, int n, sph_float dt, Params P) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i < n) sph::integrate<D>(ps[i], dt, P);
}

template <int D>
static void run(sph::Scene& s, const sph::RunConfig& cfg) {
    const int n = static_cast<int>(s.particles.size());
    const size_t bytes_p = n * sizeof(Particle), bytes_i = n * sizeof(int);
    Particle* d_ps;
    int *d_head, *d_next;
    CUDA_CHECK(cudaMalloc(&d_ps, bytes_p));
    CUDA_CHECK(cudaMalloc(&d_head, bytes_i));
    CUDA_CHECK(cudaMalloc(&d_next, bytes_i));
    CUDA_CHECK(cudaMemcpy(d_ps, s.particles.data(), bytes_p, cudaMemcpyHostToDevice));

    const int threads = 256, blocks = (n + threads - 1) / threads;
    const long steps = cfg.steps >= 0 ? cfg.steps : static_cast<long>(s.t_end / s.dt);
    const Params P = s.params;
    double compute = 0;

    if (cfg.every > 0) sph::write_vtk(cfg.out + "/" + sph::frame_filename(0), s.particles);
    for (long t = 0; t < steps; ++t) {
        auto t0 = std::chrono::steady_clock::now();
        CUDA_CHECK(cudaMemset(d_head, -1, bytes_i));
        CUDA_CHECK(cudaMemset(d_next, -1, bytes_i));
        build_cells_kernel<D><<<blocks, threads>>>(d_ps, d_head, d_next, n, P.h);
        density_pressure_kernel<D><<<blocks, threads>>>(d_ps, d_head, d_next, n, P);
        force_kernel<D><<<blocks, threads>>>(d_ps, d_head, d_next, n, P);
        integrate_kernel<D><<<blocks, threads>>>(d_ps, n, s.dt, P);
        CUDA_CHECK(cudaGetLastError());
        CUDA_CHECK(cudaDeviceSynchronize());
        compute += std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();

        if (cfg.every > 0 && (t + 1) % cfg.every == 0) {
            CUDA_CHECK(cudaMemcpy(s.particles.data(), d_ps, bytes_p, cudaMemcpyDeviceToHost));
            sph::write_vtk(cfg.out + "/" + sph::frame_filename(int(t + 1)), s.particles);
            std::cout << "step " << t + 1 << "/" << steps << "  compute " << compute << " s\n";
        }
    }
    CUDA_CHECK(cudaFree(d_ps));
    CUDA_CHECK(cudaFree(d_head));
    CUDA_CHECK(cudaFree(d_next));
    std::cout << "Simulation time (excluding I/O): " << compute << " s for " << steps << " steps\n";
}

int main(int argc, char** argv) {
    sph::RunConfig cfg = sph::parse_args(argc, argv);
    sph::Scene scene;
    if (!sph::make_scene(cfg.scene, cfg.fluid, scene)) return 1;
    if (cfg.dt > 0) scene.dt = cfg.dt;
    if (cfg.every > 0 && !sph::make_dirs(cfg.out)) return 1;

    cudaDeviceProp prop{};
    CUDA_CHECK(cudaGetDeviceProperties(&prop, 0));
    std::cout << "GPU: " << prop.name << "\n";
    sph::print_scene_summary(scene);
    if (scene.dim == 2) run<2>(scene, cfg);
    else run<3>(scene, cfg);
}
