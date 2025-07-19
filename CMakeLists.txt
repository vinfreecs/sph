cmake_minimum_required(VERSION 3.18)
project(sph LANGUAGES CUDA)

# Use C++17 standard for CUDA
set(CMAKE_CUDA_STANDARD 17)
set(CMAKE_CUDA_STANDARD_REQUIRED ON)

# Default build type
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

# Collect all CUDA source files in src/
file(GLOB_RECURSE CUDA_SOURCES CONFIGURE_DEPENDS src/*.cu)

# Add executable target
add_executable(${PROJECT_NAME} ${CUDA_SOURCES})

# Enable separate compilation (optional but recommended)
set_target_properties(${PROJECT_NAME} PROPERTIES
    CUDA_SEPARABLE_COMPILATION ON
)

# Optional: Set optimization flags explicitly for Release
# (CMake usually sets -O3 in Release mode automatically)
target_compile_options(${PROJECT_NAME} PRIVATE
  $<$<CONFIG:Release>:-O3>
)
