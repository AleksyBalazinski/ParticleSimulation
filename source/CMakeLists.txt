set(COMMON_SOURCES 
	main.cpp
	vec3.cpp
	stateRecorder.cpp
	ppMethod.cpp
	pmMethod.cpp
	simInfo.cpp
	grid.cpp
	unitConversions.cpp
	utils.cpp
	diskSampler.cpp
	diskSamplerLinear.cpp
	particle.cpp
	externalFields.cpp
	leapfrog.cpp
	chainingMesh.cpp
	p3mMethod.cpp
	greensFunctions.cpp
	demos.cpp
	FFTWAdapter.cpp
	barnesHut.cpp
)

if(USE_CUDA)
	list(APPEND COMMON_SOURCES
	printCudaVersion.cu
	cuFFTAdapter.cu
	gridGPU.cu
	leapfrogGPU.cu
	PMMethodGPU.cu
	unitConversionsGPU.cu
	)
else()
	# append CPU specific stuff
endif()

add_executable (ParticleSimulation
	${COMMON_SOURCES}
)

if(USE_CUDA)
	find_package(CUDAToolkit REQUIRED)
	set_target_properties(ParticleSimulation PROPERTIES
		CUDA_SEPARABLE_COMPILATION ON)

	target_link_libraries(ParticleSimulation PRIVATE CUDA::cufft)

endif()

set(FFTW_DIR ${PROJECT_SOURCE_DIR}/libs/fftw)



if(USE_CUDA)
	target_include_directories(ParticleSimulation
	PRIVATE ${PROJECT_SOURCE_DIR}/include
	PRIVATE ${PROJECT_SOURCE_DIR}/include_gpu
	PRIVATE ${FFTW_DIR}
	)
else()
	target_include_directories(ParticleSimulation
		PRIVATE ${PROJECT_SOURCE_DIR}/include
		PRIVATE ${FFTW_DIR}
	)
endif()

find_library(FFTW_LIB NAMES fftw3f-3 PATHS ${FFTW_DIR})
if(FFTW_LIB)
    message(STATUS "Found FFTW: ${FFTW_LIB}")
    target_link_libraries(ParticleSimulation PRIVATE ${FFTW_LIB})
else()
    message(FATAL_ERROR "FFTW library not found")
endif()

target_link_libraries(ParticleSimulation 
    PRIVATE kissfft
	PRIVATE ${FFTW_LIB}
)

add_custom_command(TARGET ParticleSimulation POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_if_different
    ${FFTW_DIR}/libfftw3f-3.dll
    $<TARGET_FILE_DIR:ParticleSimulation>)

add_library(ParticleSimulationLib)
target_sources(ParticleSimulationLib 
	PRIVATE
	vec3.cpp
	stateRecorder.cpp
	ppMethod.cpp
	barnesHut.cpp
)

target_include_directories(ParticleSimulationLib
	PUBLIC ${PROJECT_SOURCE_DIR}/include
)

target_link_libraries(ParticleSimulationLib
	PUBLIC
	kissfft
)
