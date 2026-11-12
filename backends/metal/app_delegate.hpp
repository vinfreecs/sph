#pragma once
#include <AppKit/AppKit.hpp>
#include <Metal/Metal.hpp>
#include <MetalKit/MetalKit.hpp>

class ViewDelegate;

class AppDelegate : public NS::ApplicationDelegate {
public:
    explicit AppDelegate(int fluid_particles) : fluid_(fluid_particles) {}
    ~AppDelegate() override;
    AppDelegate(const AppDelegate&) = delete;
    AppDelegate& operator=(const AppDelegate&) = delete;

    void applicationWillFinishLaunching(NS::Notification* n) override;
    void applicationDidFinishLaunching(NS::Notification* n) override;
    bool applicationShouldTerminateAfterLastWindowClosed(NS::Application*) override;

private:
    NS::Menu* createMenuBar();

    int fluid_;
    NS::Window* window_ = nullptr;
    MTK::View* view_ = nullptr;
    MTL::Device* device_ = nullptr;
    ViewDelegate* viewDelegate_ = nullptr;
};
