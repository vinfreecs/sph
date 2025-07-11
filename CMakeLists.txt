cmake_minimum_required(VERSION 3.16)
project(OpenGLApp)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_OSX_ARCHITECTURES arm64)
set(CMAKE_BUILD_TYPE Release)
# Use Homebrew's GLFW
set(GLFW_ROOT /opt/homebrew/opt/glfw)

include_directories(${GLFW_ROOT}/include)
link_directories(${GLFW_ROOT}/lib)

add_executable(OpenGLApp src/main.cpp)

target_link_libraries(OpenGLApp
    glfw
    "-framework OpenGL"
    "-framework Cocoa"
    "-framework IOKit"
    "-framework CoreVideo"
)
