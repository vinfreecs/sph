#pragma once
// Serial reference solver. Used by the CPU backend and by the OpenGL / Metal viewers.

#include <vector>

#include "physics.hpp"

namespace sph {

template <int D>
class CpuSolver {
public:
    void step(std::vector<Particle>& ps, const Params& P, sph_float dt) {
        const int n = static_cast<int>(ps.size());
        if (n == 0) return;
        head_.assign(n, -1);
        next_.assign(n, -1);
        for (int i = 0; i < n; ++i) {
            int slot = home_slot<D>(ps[i], P.h, n);
            next_[i] = head_[slot];
            head_[slot] = i;
        }
        Particle* p = ps.data();
        for (int i = 0; i < n; ++i) compute_density_pressure<D>(i, p, head_.data(), next_.data(), n, P);
        for (int i = 0; i < n; ++i) compute_forces<D>(i, p, head_.data(), next_.data(), n, P);
        for (int i = 0; i < n; ++i) integrate<D>(ps[i], dt, P);
    }

private:
    std::vector<int> head_, next_;
};

}  // namespace sph
