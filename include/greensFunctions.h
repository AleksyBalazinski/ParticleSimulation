#pragma once

#include <complex>
#include "pmConfig.h"

enum class CloudShape { S1, S2 };

std::complex<float>
GreenOptimalTSC(int kx, int ky, int kz, int dim, float a, CloudShape cs, FiniteDiffScheme fds);

std::complex<float> GreenDiscreteLaplacian(int kx, int ky, int kz, int dim);