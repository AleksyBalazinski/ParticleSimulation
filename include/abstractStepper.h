#pragma once

#include <vector>
template <typename T>
class AbstractStepper {
 public:
  virtual ~AbstractStepper() {}
  virtual void doStep(std::vector<T>& x, float dt) = 0;
};
