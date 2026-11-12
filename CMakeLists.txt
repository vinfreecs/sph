cmake_minimum_required(VERSION 3.18)
project(sph LANGUAGES CXX)

option(SPH_CPU    "Serial CPU reference solver (sph_cpu)"     ON)
option(SPH_CUDA   "CUDA solver (sph_cuda)"                    OFF)
option(SPH_OPENGL "2D OpenGL viewer, needs GLFW (sph_opengl)" OFF)
option(SPH_METAL  "2D Metal viewer, macOS only (sph_metal)"   OFF)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

# Header-only physics shared by every backend
add_library(sph_core INTERFACE)
target_include_directories(sph_core INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/core/include)

if(SPH_CPU)
  add_executable(sph_cpu backends/cpu/main.cpp)
  target_link_libraries(sph_cpu PRIVATE sph_core)
endif()

if(SPH_CUDA)
  # A100 = sm_80, A40 = sm_86 (NHR@FAU Alex)
  if(NOT DEFINED CMAKE_CUDA_ARCHITECTURES)
    set(CMAKE_CUDA_ARCHITECTURES 80 86)
  endif()
  enable_language(CUDA)
  set(CMAKE_CUDA_STANDARD 17)
  set(CMAKE_CUDA_STANDARD_REQUIRED ON)
  add_executable(sph_cuda backends/cuda/main.cu)
  target_link_libraries(sph_cuda PRIVATE sph_core)
endif()

if(SPH_OPENGL)
  find_package(glfw3 3.3 REQUIRED)
  find_package(OpenGL REQUIRED)
  add_library(glad STATIC third_party/glad/src/glad.c)
  set_target_properties(glad PROPERTIES LINKER_LANGUAGE C)
  enable_language(C)
  target_include_directories(glad PUBLIC third_party/glad/include)
  add_executable(sph_opengl backends/opengl/main.cpp)
  target_link_libraries(sph_opengl PRIVATE sph_core glad glfw OpenGL::GL)
endif()

if(SPH_METAL)
  if(NOT APPLE)
    message(FATAL_ERROR "SPH_METAL requires macOS")
  endif()
  enable_language(OBJCXX)
  add_executable(sph_metal
    backends/metal/main.cpp
    backends/metal/app_delegate.cpp
    backends/metal/view_delegate.cpp
    backends/metal/renderer.cpp
    backends/metal/metal_cpp_impl.mm)
  target_include_directories(sph_metal SYSTEM PRIVATE third_party/metal-cpp)
  # metal-cpp's implementation unit must be compiled without ARC
  set_source_files_properties(backends/metal/metal_cpp_impl.mm PROPERTIES COMPILE_OPTIONS -fno-objc-arc)
  target_link_libraries(sph_metal PRIVATE sph_core
    "-framework Foundation" "-framework AppKit" "-framework Metal"
    "-framework MetalKit" "-framework QuartzCore")
endif()
