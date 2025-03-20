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

float TSCFourier(std::array<float, 3> k) {
  float prod = 1;
  for (int i = 0; i < 3; i++) {
    prod *= std::powf(sinc(k[i] / 2), 3);
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

std::array<std::complex<float>, 3> RFourier(std::array<float, 3> k, float a, CloudShape cs) {
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

std::array<std::complex<float>, 3> D2Fourier(std::array<float, 3> k) {
  std::array<std::complex<float>, 3> D;
  std::complex<float> I(0, 1);
  for (int i = 0; i < 3; i++) {
    D[i] = I * std::sinf(k[i]);
  }

  return D;
}

std::complex<float> dotProduct(std::array<std::complex<float>, 3> a,
                               std::array<std::complex<float>, 3> b) {
  std::complex<float> dot(0, 0);
  for (int i = 0; i < 3; i++) {
    std::complex<float> biConj(b[i].real(), -1 * b[i].imag());
    dot += a[i] * biConj;
  }

  return dot;
}

float TSCAliasSum(std::array<float, 3> k) {
  float sum = 1;
  for (int i = 0; i < 3; i++) {
    sum *= (1 - std::powf(std::sinf(k[i] / 2), 2) + 2.0f / 15 * std::powf(std::sinf(k[i] / 2), 4));
  }
  return sum;
}

std::complex<float> GreenOptimalTSC(int kx,
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

  float denomSum = TSCAliasSum(k);

  std::array<std::complex<float>, 3> D;
  if (fds == FiniteDiffScheme::TWO_POINT) {
    D = D2Fourier(k);
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
        float uSquared = std::powf(TSCFourier(kn), 2);

        std::array<std::complex<float>, 3> R = RFourier(kn, a, cs);
        numDotRight[0] += uSquared * R[0];
        numDotRight[1] += uSquared * R[1];
        numDotRight[2] += uSquared * R[2];
      }
    }
  }

  std::complex<float> numerator = dotProduct(D, numDotRight);
  float denominator = DNormSquared * denomSum;

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
