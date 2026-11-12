#pragma once
#include <MetalKit/MetalKit.hpp>

#include "renderer.hpp"
#include "sph/app.hpp"
#include "sph/cpu_solver.hpp"

// Advances the simulation one step per frame and draws it.
class ViewDelegate : public MTK::ViewDelegate {
public:
    ViewDelegate(MTL::Device* device, sph::Scene scene);
    ~ViewDelegate() override;
    void drawInMTKView(MTK::View* view) override;

private:
    Renderer* renderer_ = nullptr;
    sph::Scene scene_;
    sph::CpuSolver<2> solver_;
};
