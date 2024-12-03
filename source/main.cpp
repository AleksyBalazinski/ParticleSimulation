#include <chrono>
#include <cmath>
#include <iostream>
#include <string>
#include "pmMethod.h"
#include "ppMethod.h"
#include "utils.h"

int main() {
  auto posRange = std::make_pair(Vec3(5.0, 5.0, 5.0), Vec3(20.0, 30.0, 30.0));
  auto vRange = std::make_pair(Vec3(-1.0, -1.0, -1.0), Vec3(1.0, 1.0, 1.0));
  auto massRange = std::make_pair(20, 30);

  const int N = 3;
  auto state = randomInitialState(N, posRange, vRange);
  auto masses = randomMasses(N, massRange);

  PMMethod pm(state, masses, 32, 1, 0.001);
  pm.run(3.0);
  // ppMethod(state, masses, 3.0);

  // ppMethod(state, masses, 3);

  // auto posRange = std::make_pair(Vec3(-10.0, -10.0, -10.0), Vec3(10.0, 10.0, 10.0));
  // auto vRange = std::make_pair(Vec3(-2.0, -2.0, -2.0), Vec3(2.0, 2.0, 0.0));
  // auto massRange = std::make_pair(20, 30);

  // const int N = 3;
  // auto state = randomInitialState(N, posRange, vRange);
  // auto masses = randomMasses(N, massRange);

  // std::string saveDirPath = ppMethod(state, masses, 10.0);
  // std::cout << "saving at " << saveDirPath << '\n';

  // std::cout << "2D FFT\n";
  // float input[] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18};
  // int nfft = sizeof(input) / sizeof(float);  // nfft = 8
  // kiss_fft_cpx* cin = new kiss_fft_cpx[nfft];
  // kiss_fft_cpx* cout = new kiss_fft_cpx[nfft];
  // cin = new kiss_fft_cpx[nfft];
  // cout = new kiss_fft_cpx[nfft];

  // memset(cin, 0, nfft * sizeof(kiss_fft_cpx));
  // memset(cout, 0, nfft * sizeof(kiss_fft_cpx));

  // for (int i = 0; i < nfft; i++) {
  //   cin[i].r = input[i];
  // }

  // bool is_inverse_fft = false;
  // int dims[] = {3, 3, 2};
  // kiss_fftnd_cfg cfg_ndf = kiss_fftnd_alloc(dims, 3, is_inverse_fft, 0, 0);

  // kiss_fftnd(cfg_ndf, cin, cout);

  // printf("\nForward Transform:\n");
  // for (int i = 0; i < nfft; ++i) {
  //   printf("#%d  %f %fj\n", i, cout[i].r, cout[i].i);
  // }

  // is_inverse_fft = true;
  // kiss_fftnd_cfg cfg_ndf_i = kiss_fftnd_alloc(dims, 3, is_inverse_fft, 0, 0);

  // kiss_fftnd(cfg_ndf_i, cout, cin);

  // // original input data
  // printf("\nInverse Transform:\n");
  // for (int i = 0; i < nfft; ++i) {
  //   printf("#%d  %f\n", i, cin[i].r / nfft);  // div by N to scale data back to the original
  //   range
  // }

  // kiss_fft_free(cfg_ndf);
  // kiss_fft_free(cfg_ndf_i);
  // delete[] cin;
  // delete[] cout;
}