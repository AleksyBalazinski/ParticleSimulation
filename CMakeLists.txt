cmake_minimum_required(VERSION 3.13)

project(ParticleSimulation LANGUAGES CXX)

# set the C++ standard
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

add_subdirectory(libs/kissfft)
add_subdirectory(source)

enable_testing()

add_subdirectory(test)
