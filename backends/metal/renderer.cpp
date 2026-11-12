#include "renderer.hpp"

#include <simd/simd.h>

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstring>
#include <limits>

namespace {

const char* kShaderSrc = R"(
#include <metal_stdlib>
using namespace metal;

struct Instance {
    float2 pos;      // world position, the scene lives in roughly [-1,1]
    float  radius;
    float3 color;
    float  fixedFlag;
};
struct Params { float aspect; };
struct VSOut {
    float4 position [[position]];
    float2 local;
    float3 color;
};

vertex VSOut v_main(uint vid [[vertex_id]], uint iid [[instance_id]],
                    constant float2* quad [[buffer(0)]],
                    device const Instance* inst [[buffer(1)]],
                    constant Params& P [[buffer(2)]]) {
    Instance s = inst[iid];
    float2 q = quad[vid];
    // keep circles round and the world square visible for any window aspect
    float2 scale = P.aspect >= 1.0 ? float2(1.0 / P.aspect, 1.0) : float2(1.0, P.aspect);
    VSOut o;
    o.position = float4((s.pos + q * s.radius) * scale, 0.0, 1.0);
    o.local = q;
    o.color = s.fixedFlag > 0.5 ? float3(1.0) : s.color;
    return o;
}

fragment half4 f_main(VSOut in [[stage_in]]) {
    if (length(in.local) > 1.0) discard_fragment();
    return half4(half3(in.color), 1.0h);
}
)";

struct InstanceCPU {  // must match `Instance` above (float3 is 16-byte aligned)
    simd::float2 pos;
    float radius;
    simd::float3 color;
    float fixedFlag;
};

void speedToColor(float t, float& r, float& g, float& b) {  // blue -> cyan -> green -> yellow -> red
    t = std::clamp(t, 0.0f, 1.0f);
    if (t < 0.25f) { r = 0; g = t * 4; b = 1; }
    else if (t < 0.5f) { r = 0; g = 1; b = 1 - (t - 0.25f) * 4; }
    else if (t < 0.75f) { r = (t - 0.5f) * 4; g = 1; b = 0; }
    else { r = 1; g = 1 - (t - 0.75f) * 4; b = 0; }
}

}  // namespace

Renderer::Renderer(MTL::Device* device) : device_(device->retain()) {
    queue_ = device_->newCommandQueue();
    buildPipeline();
}

Renderer::~Renderer() {
    if (instanceBuf_) instanceBuf_->release();
    if (quadVB_) quadVB_->release();
    if (pso_) pso_->release();
    if (queue_) queue_->release();
    if (device_) device_->release();
}

void Renderer::buildPipeline() {
    using NS::StringEncoding::UTF8StringEncoding;
    NS::Error* err = nullptr;
    MTL::Library* lib = device_->newLibrary(NS::String::string(kShaderSrc, UTF8StringEncoding), nullptr, &err);
    assert(lib && "particle shader failed to compile");
    MTL::Function* vfn = lib->newFunction(NS::String::string("v_main", UTF8StringEncoding));
    MTL::Function* ffn = lib->newFunction(NS::String::string("f_main", UTF8StringEncoding));

    MTL::RenderPipelineDescriptor* desc = MTL::RenderPipelineDescriptor::alloc()->init();
    desc->setVertexFunction(vfn);
    desc->setFragmentFunction(ffn);
    desc->colorAttachments()->object(0)->setPixelFormat(MTL::PixelFormatBGRA8Unorm_sRGB);
    pso_ = device_->newRenderPipelineState(desc, &err);
    assert(pso_ && "pipeline creation failed");
    desc->release();
    vfn->release();
    ffn->release();
    lib->release();

    const simd::float2 quad[6] = {{-1, -1}, {1, -1}, {-1, 1}, {1, -1}, {1, 1}, {-1, 1}};
    quadVB_ = device_->newBuffer(quad, sizeof quad, MTL::ResourceStorageModeShared);
}

void Renderer::drawParticles(MTK::View* view, const std::vector<sph::Particle>& particles) {
    const size_t count = particles.size();
    if (count == 0) return;

    float vmin = std::numeric_limits<float>::max(), vmax = 0;
    for (const auto& p : particles) {
        float v = std::sqrt(p.vel_x * p.vel_x + p.vel_y * p.vel_y);
        vmin = std::min(vmin, v);
        vmax = std::max(vmax, v);
    }
    const float inv_range = 1.0f / std::max(1e-5f, vmax - vmin);

    const size_t bytes = count * sizeof(InstanceCPU);
    if (!instanceBuf_ || instanceBuf_->length() < bytes) {
        if (instanceBuf_) instanceBuf_->release();
        instanceBuf_ = device_->newBuffer(bytes, MTL::ResourceStorageModeShared);
    }
    auto* dst = static_cast<InstanceCPU*>(instanceBuf_->contents());
    for (size_t i = 0; i < count; ++i) {
        const auto& p = particles[i];
        float r = 1, g = 1, b = 1;
        if (!p.fixed) speedToColor((std::sqrt(p.vel_x * p.vel_x + p.vel_y * p.vel_y) - vmin) * inv_range, r, g, b);
        dst[i] = {{p.pos_x, p.pos_y}, p.radius, {r, g, b}, p.fixed ? 1.0f : 0.0f};
    }

    CGSize size = view->drawableSize();
    float aspect = size.height > 0 ? float(size.width / size.height) : 1.0f;

    NS::AutoreleasePool* pool = NS::AutoreleasePool::alloc()->init();
    MTL::RenderPassDescriptor* rpd = view->currentRenderPassDescriptor();
    if (rpd) {
        MTL::CommandBuffer* cmd = queue_->commandBuffer();
        MTL::RenderCommandEncoder* enc = cmd->renderCommandEncoder(rpd);
        enc->setRenderPipelineState(pso_);
        enc->setVertexBuffer(quadVB_, 0, 0);
        enc->setVertexBuffer(instanceBuf_, 0, 1);
        enc->setVertexBytes(&aspect, sizeof aspect, 2);
        enc->drawPrimitives(MTL::PrimitiveTypeTriangle, NS::UInteger(0), NS::UInteger(6), NS::UInteger(count));
        enc->endEncoding();
        cmd->presentDrawable(view->currentDrawable());
        cmd->commit();
    }
    pool->release();
}
