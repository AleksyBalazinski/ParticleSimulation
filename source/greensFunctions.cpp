#include "greensFunctions.h"

#include <array>
#include <numbers>

constexpr float pi = std::numbers::pi_v<float>;

float sinc(float x) {
  if (x == 0) {
    return 1;
  }
  return std::sinf(x) / x;
}

float TSCFourier(const std::array<float, 3>& k) {
  float prod = 1;
  for (int i = 0; i < 3; i++) {
    prod *= std::powf(sinc(k[i] / 2), 3);
  }

  return prod;
}

float CICFourier(const std::array<float, 3>& k) {
  float prod = 1;
  for (int i = 0; i < 3; i++) {
    prod *= std::powf(sinc(k[i] / 2), 2);
  }

  return prod;
}

float NGPFourier(const std::array<float, 3>& k) {
  float prod = 1;
  for (int i = 0; i < 3; i++) {
    prod *= sinc(k[i] / 2);
  }

  return prod;
}

float S1Fourier(float k, float a) {
  float u = k * a / 2;
  return -3 / std::powf(u, 3) * (u * std::cosf(u) - std::sinf(u));
}

float S2Fourier(float k, float a) {
  float u = k * a / 2;
  return 12 / std::powf(u, 4) * (2 - 2 * std::cosf(u) - u * std::sinf(u));
}

std::array<std::complex<float>, 3> RFourier(const std::array<float, 3>& k, float a, CloudShape cs) {
  std::array<std::complex<float>, 3> R;
  std::complex<float> I(0, 1);
  float kLength = std::sqrtf(k[0] * k[0] + k[1] * k[1] + k[2] * k[2]);
  float s;
  if (cs == CloudShape::S1) {
    s = S1Fourier(kLength, a);
  } else if (cs == CloudShape::S2) {
    s = S2Fourier(kLength, a);
  }
  float sSquared = std::powf(s, 2);
  for (int i = 0; i < 3; i++) {
    R[i] = -I * k[i] * sSquared / (kLength * kLength);
  }

  return R;
}

std::array<std::complex<float>, 3> D2Fourier(const std::array<float, 3>& k) {
  std::array<std::complex<float>, 3> D;
  std::complex<float> I(0, 1);
  for (int i = 0; i < 3; i++) {
    D[i] = I * std::sinf(k[i]);
  }

  return D;
}

std::array<std::complex<float>, 3> D4Fourier(const std::array<float, 3>& k) {
  std::array<std::complex<float>, 3> D;
  std::complex<float> I(0, 1);
  float alpha = 4.0f / 3;
  for (int i = 0; i < 3; i++) {
    D[i] = I * alpha * std::sinf(k[i]) + I * (1 - alpha) * std::sinf(2 * k[i]) / 2.0f;
  }

  return D;
}

std::complex<float> dotProduct(const std::array<std::complex<float>, 3>& a,
                               const std::array<std::complex<float>, 3>& b) {
  std::complex<float> dot(0, 0);
  for (int i = 0; i < 3; i++) {
    std::complex<float> biConj(b[i].real(), -1 * b[i].imag());
    dot += a[i] * biConj;
  }

  return dot;
}

float TSCAliasSum(const std::array<float, 3>& k) {
  float sum = 1;
  for (int i = 0; i < 3; i++) {
    sum *= (1 - std::powf(std::sinf(k[i] / 2), 2) + 2.0f / 15 * std::powf(std::sinf(k[i] / 2), 4));
  }
  return sum;
}

float CICAliasSum(const std::array<float, 3>& k) {
  float sum = 1;
  for (int i = 0; i < 3; i++) {
    sum *= (1 + 2 * std::powf(std::cosf(k[i] / 2), 2));
  }
  return (1.0f / 27) * sum;
}

float NGPAliasSum(const std::array<float, 3>& k) {
  return 1;
}

std::complex<float> GreenOptimal(InterpolationScheme is,
                                 int kx,
                                 int ky,
                                 int kz,
                                 std::tuple<int, int, int> dims,
                                 float a,
                                 CloudShape cs,
                                 FiniteDiffScheme fds) {
  if (kx == 0 && ky == 0 && kz == 0) {
    return 0;
  }

  std::array<float, 3> k = {2 * pi * float(kx) / std::get<0>(dims),
                            2 * pi * float(ky) / std::get<1>(dims),
                            2 * pi * float(kz) / std::get<2>(dims)};

  float denomSum;
  if (is == InterpolationScheme::TSC) {
    denomSum = TSCAliasSum(k);
  } else if (is == InterpolationScheme::CIC) {
    denomSum = CICAliasSum(k);
  } else if (is == InterpolationScheme::NGP) {
    denomSum = NGPAliasSum(k);
  } else {
    throw std::invalid_argument("Not implemented");
  }

  std::array<std::complex<float>, 3> D;
  if (fds == FiniteDiffScheme::TWO_POINT) {
    D = D2Fourier(k);
  } else if (fds == FiniteDiffScheme::FOUR_POINT) {
    D = D4Fourier(k);
  } else {
    throw std::invalid_argument("not implemented");
  }

  float DNormSquared = dotProduct(D, D).real();

  int limit = 2;
  std::array<std::complex<float>, 3> numDotRight = {0, 0, 0};
  for (int n1 = -limit; n1 <= limit; n1++) {
    for (int n2 = -limit; n2 <= limit; n2++) {
      for (int n3 = -limit; n3 <= limit; n3++) {
        std::array<float, 3> kn = {k[0] + 2 * pi * n1, k[1] + 2 * pi * n2, k[2] + 2 * pi * n3};
        float uSquared;
        if (is == InterpolationScheme::TSC) {
          uSquared = std::powf(TSCFourier(kn), 2);
        } else if (is == InterpolationScheme::CIC) {
          uSquared = std::powf(CICFourier(kn), 2);
        } else if (is == InterpolationScheme::NGP) {
          uSquared = std::powf(NGPFourier(kn), 2);
        } else {
          throw std::invalid_argument("Not implemented");
        }

        std::array<std::complex<float>, 3> R = RFourier(kn, a, cs);
        numDotRight[0] += uSquared * R[0];
        numDotRight[1] += uSquared * R[1];
        numDotRight[2] += uSquared * R[2];
      }
    }
  }

  std::complex<float> numerator = dotProduct(D, numDotRight);
  float denominator = DNormSquared * denomSum * denomSum;

  return numerator / denominator;
}

std::complex<float> GreenDiscreteLaplacian(int kx, int ky, int kz, std::tuple<int, int, int> dims) {
  if (kx == 0 && ky == 0 && kz == 0) {
    return 0;
  }

  auto sx = std::sinf(pi * kx / std::get<0>(dims));
  auto sy = std::sinf(pi * ky / std::get<1>(dims));
  auto sz = std::sinf(pi * kz / std::get<2>(dims));
  return -0.25f / (sx * sx + sy * sy + sz * sz);
}

std::complex<float> GreenPoorMan(int i, int j, int k, std::tuple<int, int, int> dims) {
  if (i == 0 && j == 0 && k == 0) {
    return 0.0f;
  }

  auto pi = std::numbers::pi_v<float>;

  auto [Nx, Ny, Nz] = dims;
  int ki = (i <= Nx / 2) ? i : i - Nx;
  int kj = (j <= Ny / 2) ? j : j - Ny;
  int kk = (k <= Nz / 2) ? k : k - Nz;

  float kx = 2 * pi * ki / Nx;
  float ky = 2 * pi * kj / Ny;
  float kz = 2 * pi * kk / Nz;

  float k2 = kx * kx + ky * ky + kz * kz;

  return -1.0f / k2;
}