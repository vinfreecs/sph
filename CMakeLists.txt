# cmake_minimum_required(VERSION 3.18)
# project(sph LANGUAGES CUDA)

# # Set the CUDA standard
# set(CMAKE_CUDA_STANDARD 17)
# set(CMAKE_CUDA_STANDARD_REQUIRED ON)

# # Set Release build by default
# set(CMAKE_BUILD_TYPE Release)

# # Get all .cu files from src/
# file(GLOB_RECURSE CUDA_SOURCES CONFIGURE_DEPENDS src/*.cu)

# # Create the CUDA executable
# add_executable(${PROJECT_NAME} ${CUDA_SOURCES})

# # Optional: Enable separate compilation for device code
# set_target_properties(${PROJECT_NAME} PROPERTIES
#     CUDA_SEPARABLE_COMPILATION ON
# )


cmake_minimum_required(VERSION 3.12)
project(sph LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CMAKE_BUILD_TYPE Release)
# Add all .cpp files in src/ as source files
file(GLOB_RECURSE SOURCES "src/*.cpp")
# Define the executable
add_executable(${PROJECT_NAME} ${SOURCES})

