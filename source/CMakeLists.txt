add_executable (ParticleSimulation
	main.cpp
	vec3.cpp
	stateRecorder.cpp
	ppMethod.cpp
	pmMethod.cpp
	simInfo.cpp
	grid.cpp
	unitConversions.cpp
)

target_include_directories(ParticleSimulation
	PRIVATE ${PROJECT_SOURCE_DIR}/include
)

target_link_libraries(ParticleSimulation 
    PRIVATE kissfft
    PRIVATE glfw
    PRIVATE libglew_static
    PRIVATE glm
    PRIVATE pocketfft
)

add_library(ParticleSimulationLib)
target_sources(ParticleSimulationLib 
	PRIVATE
	vec3.cpp
	stateRecorder.cpp
	ppMethod.cpp
)

target_include_directories(ParticleSimulationLib
	PUBLIC ${PROJECT_SOURCE_DIR}/include
)

target_link_libraries(ParticleSimulationLib
	PUBLIC
	pocketfft
	kissfft
)
