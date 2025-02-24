NVCC := nvcc
CXXFLAGS := -std=c++17 -dc
LDFLAGS := -lcufft

SRC := conversions.cu disk_sampler_linear.cu ext_field_bind.cu main.cu pm.cu state_recorder.cu utils.cu vec3.cu
OBJ := bin/conversions.obj bin/disk_sampler_linear.obj bin/ext_field_bind.obj bin/main.obj bin/pm.obj bin/state_recorder.obj bin/utils.obj bin/vec3.obj
BIN_DIR := bin
EXECUTABLE := $(BIN_DIR)/pm_gpu


$(EXECUTABLE): $(BIN_DIR) $(OBJ)
	$(NVCC) $(OBJ) -o $@ $(LDFLAGS)

$(BIN_DIR):
	powershell -Command "if (!(Test-Path '$(BIN_DIR)')) { New-Item -ItemType Directory -Path '$(BIN_DIR)' }"

bin/%.obj: %.cu
	$(NVCC) $(CXXFLAGS) -c $< -o $@

