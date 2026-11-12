#include "view_delegate.hpp"

#include <utility>

ViewDelegate::ViewDelegate(MTL::Device* device, sph::Scene scene)
    : MTK::ViewDelegate(), renderer_(new Renderer(device)), scene_(std::move(scene)) {}

ViewDelegate::~ViewDelegate() { delete renderer_; }

void ViewDelegate::drawInMTKView(MTK::View* view) {
    solver_.step(scene_.particles, scene_.params, scene_.dt);
    renderer_->drawParticles(view, scene_.particles);
}
