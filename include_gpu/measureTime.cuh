#pragma once

#include <chrono>

#define declareDeviceTimer(event)  \
  cudaEvent_t cuda##event##Start_; \
  cudaEvent_t cuda##event##Stop_;  \
  float cuda##event##TimeMs_ = 0;

#define declareHostTimer(event) float host##event##TimeMs_ = 0;

#define measureDeviceTime(event, ...)     \
  do {                                    \
    cudaEventRecord(cuda##event##Start_); \
    __VA_ARGS__;                          \
    cudaEventRecord(cuda##event##Stop_);  \
  } while (0)

#define accDeviceTime(event)                                              \
  do {                                                                    \
    cudaEventSynchronize(cuda##event##Stop_);                             \
    float acc_ = 0;                                                       \
    cudaEventElapsedTime(&acc_, cuda##event##Start_, cuda##event##Stop_); \
    cuda##event##TimeMs_ += acc_;                                         \
  } while (0)

#define createCudaEvents(event)          \
  cudaEventCreate(&cuda##event##Start_); \
  cudaEventCreate(&cuda##event##Stop_);

#define destroyCudaEvents(event)         \
  cudaEventDestroy(cuda##event##Start_); \
  cudaEventDestroy(cuda##event##Stop_);

#define measureHostTime(event, ...)                               \
  do {                                                            \
    auto event##Start_ = std::chrono::steady_clock::now();        \
    __VA_ARGS__;                                                  \
    auto event##Stop_ = std::chrono::steady_clock::now();         \
    host##event##Time##Ms_ += toMs(event##Stop_ - event##Start_); \
  } while (0)

#define printHostTime(event)                                          \
  do {                                                                \
    std::cout << #event << ": " << host##event##Time##Ms_ << " ms\n"; \
  } while (0)

#define printDeviceTime(event)                                        \
  do {                                                                \
    std::cout << #event << ": " << cuda##event##Time##Ms_ << " ms\n"; \
  } while (0)

inline float toMs(const std::chrono::nanoseconds delta) {
  return delta.count() * 1e-6f;
}