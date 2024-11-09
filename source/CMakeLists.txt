add_executable (ParticleSimulation
	main.cpp
	abstractStepper.h
	RK4Stepper.h
	vec3.h
	vec3.cpp
	utils.h
	stateRecorder.h
	stateRecorder.cpp
	ppMethod.h 
	ppMethod.cpp
	)

add_library(ParticleSimulationLib)
target_sources(ParticleSimulationLib 
PRIVATE
	vec3.cpp
	stateRecorder.cpp
	ppMethod.cpp
)

target_include_directories(ParticleSimulationLib
PUBLIC 
	${CMAKE_CURRENT_LIST_DIR}
)
