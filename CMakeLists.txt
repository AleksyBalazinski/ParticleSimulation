cmake_minimum_required(VERSION 3.13)

project(ParticleSimulation LANGUAGES CXX)

# set the C++ standard
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

add_subdirectory(libs/kissfft EXCLUDE_FROM_ALL)
add_subdirectory(libs/glfw EXCLUDE_FROM_ALL)
add_subdirectory(libs/glew EXCLUDE_FROM_ALL)
add_subdirectory(libs/glm EXCLUDE_FROM_ALL)
add_subdirectory(libs/pocketfft EXCLUDE_FROM_ALL)
add_subdirectory(source)

target_link_libraries(ParticleSimulation 
    PRIVATE kissfft
    PRIVATE glfw
    PRIVATE libglew_static
    PRIVATE glm
    PRIVATE pocketfft
)

enable_testing()

add_subdirectory(test)
