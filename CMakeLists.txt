cmake_minimum_required(VERSION 3.10)

project(ParticleSimulation LANGUAGES CXX)

option(USE_CUDA "Enable CUDA support" OFF)

if(USE_CUDA)
    message(STATUS "CUDA enabled")
    enable_language("CUDA")
    set(CMAKE_CUDA_STANDARD 20)
    set(CMAKE_CUDA_STANDARD_REQUIRED ON)
    add_definitions(-DCUDA)
endif()

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

add_subdirectory(libs/kissfft EXCLUDE_FROM_ALL)
add_subdirectory(libs/glfw EXCLUDE_FROM_ALL)
add_subdirectory(libs/glew EXCLUDE_FROM_ALL)
add_subdirectory(libs/glm EXCLUDE_FROM_ALL)
add_subdirectory(source)

enable_testing()

add_subdirectory(test)
