// 2D Metal viewer (metal-cpp): same "cup" scene and CPU solver as the OpenGL viewer.

#include <AppKit/AppKit.hpp>

#include <cstdlib>

#include "app_delegate.hpp"

int main(int argc, char** argv) {
    int fluid = argc > 1 ? std::atoi(argv[1]) : 1000;
    NS::AutoreleasePool* pool = NS::AutoreleasePool::alloc()->init();
    AppDelegate delegate(fluid);
    NS::Application* app = NS::Application::sharedApplication();
    app->setDelegate(&delegate);
    app->run();
    pool->release();
    return 0;
}
