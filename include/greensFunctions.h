#pragma once

#include <complex>
#include "pmConfig.h"

std::complex<float> GreenS1OptimalTSC(int kx,
                                      int ky,
                                      int kz,
                                      int dim,
                                      float a,
                                      FiniteDiffScheme fds);

std::complex<float> GreenDiscreteLaplacian(int kx, int ky, int kz, int dim);