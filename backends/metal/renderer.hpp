#pragma once
#include <Metal/Metal.hpp>
#include <MetalKit/MetalKit.hpp>
#include <vector>

#include "sph/particle.hpp"

// Draws particles as instanced, analytically-masked quads coloured by speed.
class Renderer {
public:
    explicit Renderer(MTL::Device* device);
    ~Renderer();
    void drawParticles(MTK::View* view, const std::vector<sph::Particle>& particles);

private:
    void buildPipeline();

    MTL::Device* device_ = nullptr;
    MTL::CommandQueue* queue_ = nullptr;
    MTL::RenderPipelineState* pso_ = nullptr;
    MTL::Buffer* quadVB_ = nullptr;
    MTL::Buffer* instanceBuf_ = nullptr;
};
