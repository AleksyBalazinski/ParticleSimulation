#pragma once

#include <complex>
#include <tuple>
#include "pmConfig.h"

enum class CloudShape { S1, S2 };

std::complex<float> GreenOptimalTSC(int kx,
                                    int ky,
                                    int kz,
                                    std::tuple<int, int, int> dims,
                                    float a,
                                    CloudShape cs,
                                    FiniteDiffScheme fds);

std::complex<float> GreenDiscreteLaplacian(int kx, int ky, int kz, std::tuple<int, int, int> dims);