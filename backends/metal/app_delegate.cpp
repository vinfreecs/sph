#include "app_delegate.hpp"

#include <cstdlib>
#include <utility>

#include "view_delegate.hpp"

using NS::StringEncoding::UTF8StringEncoding;

AppDelegate::~AppDelegate() {
    if (view_) view_->release();
    if (window_) window_->release();
    if (device_) device_->release();
    delete viewDelegate_;
}

NS::Menu* AppDelegate::createMenuBar() {
    NS::Menu* mainMenu = NS::Menu::alloc()->init();
    NS::MenuItem* appItem = NS::MenuItem::alloc()->init();
    NS::Menu* appMenu = NS::Menu::alloc()->init(NS::String::string("SPH", UTF8StringEncoding));
    NS::String* quitName = NS::String::string("Quit ", UTF8StringEncoding)
                               ->stringByAppendingString(NS::RunningApplication::currentApplication()->localizedName());
    SEL quit = NS::MenuItem::registerActionCallback(
        "appQuit", [](void*, SEL, const NS::Object* sender) { NS::Application::sharedApplication()->terminate(sender); });
    NS::MenuItem* quitItem = appMenu->addItem(quitName, quit, NS::String::string("q", UTF8StringEncoding));
    quitItem->setKeyEquivalentModifierMask(NS::EventModifierFlagCommand);
    appItem->setSubmenu(appMenu);
    mainMenu->addItem(appItem);
    appItem->release();
    appMenu->release();
    return mainMenu->autorelease();
}

void AppDelegate::applicationWillFinishLaunching(NS::Notification* n) {
    auto* app = reinterpret_cast<NS::Application*>(n->object());
    app->setMainMenu(createMenuBar());
    app->setActivationPolicy(NS::ActivationPolicy::ActivationPolicyRegular);
}

void AppDelegate::applicationDidFinishLaunching(NS::Notification* n) {
    CGRect frame = {{100, 100}, {1200, 800}};
    window_ = NS::Window::alloc()->init(frame, NS::WindowStyleMaskClosable | NS::WindowStyleMaskTitled,
                                        NS::BackingStoreBuffered, false);
    device_ = MTL::CreateSystemDefaultDevice();
    view_ = MTK::View::alloc()->init(frame, device_);
    view_->setColorPixelFormat(MTL::PixelFormatBGRA8Unorm_sRGB);
    view_->setClearColor(MTL::ClearColor::Make(0.1, 0.1, 0.1, 1.0));

    sph::Scene scene;
    if (!sph::make_scene("cup", fluid_, scene)) std::exit(1);
    sph::print_scene_summary(scene);
    viewDelegate_ = new ViewDelegate(device_, std::move(scene));
    view_->setDelegate(viewDelegate_);

    window_->setContentView(view_);
    window_->setTitle(NS::String::string("SPH (Metal)", UTF8StringEncoding));
    window_->makeKeyAndOrderFront(nullptr);
    reinterpret_cast<NS::Application*>(n->object())->activateIgnoringOtherApps(true);
}

bool AppDelegate::applicationShouldTerminateAfterLastWindowClosed(NS::Application*) { return true; }
