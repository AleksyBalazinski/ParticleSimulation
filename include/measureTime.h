#pragma once

#include <chrono>
#include <iostream>

#define declareTimeAcc(event) float event##TimeMs_ = 0;

#define measureTime(event, ...)                                       \
  do {                                                                \
    auto event##Start_ = std::chrono::steady_clock::now();            \
    __VA_ARGS__;                                                      \
    auto event##Stop_ = std::chrono::steady_clock::now();             \
    event##TimeMs_ += (event##Stop_ - event##Start_).count() * 1e-6f; \
  } while (0)

#define printTime(event)                                        \
  do {                                                          \
    std::cout << #event << ": " << event##Time##Ms_ << " ms\n"; \
  } while (0)