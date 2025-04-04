NVCC := nvcc
CXXFLAGS := -std=c++17 -dc -Wno-deprecated-gpu-targets
LDFLAGS := -lcufft

BIN_DIR := bin
EXECUTABLE := $(BIN_DIR)/pm_gpu

SRC := conversions.cu disk_sampler_linear.cu main.cu pm.cu state_recorder.cu utils.cu vec3.cu external_fields.cu
_OBJ := conversions.obj disk_sampler_linear.obj main.obj pm.obj state_recorder.obj utils.obj vec3.obj external_fields.obj
OBJ := $(patsubst %, $(BIN_DIR)/%, $(_OBJ))
DEPS := settings.cuh conversions.cuh disk_sampler_linear.cuh external_fields.cuh helper_macros.h pm.cuh state_recorder.cuh utils.cuh vec3.cuh

$(EXECUTABLE): $(BIN_DIR) $(OBJ)
	$(NVCC) $(OBJ) -o $@ $(LDFLAGS)

$(BIN_DIR):
	powershell -Command "if (!(Test-Path '$(BIN_DIR)')) { New-Item -ItemType Directory -Path '$(BIN_DIR)' }"

bin/%.obj: %.cu $(DEPS)
	$(NVCC) $(CXXFLAGS) -c $< -o $@

