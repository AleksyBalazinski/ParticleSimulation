#include <chrono>
#include <cmath>
#include <iostream>
#include <string>
#include "kissFFTAdapter.h"
#include "pmMethod.h"
#include "ppMethod.h"
#include "utils.h"

#include <GL/glew.h>
#include <GLFW/glfw3.h>

void framebufferSizeCallback(GLFWwindow* window, int width, int height);

int main() {
  // const int N = 2;
  // std::vector<Vec3> state = {Vec3(10, 10, 10), Vec3(25, 15, 20), Vec3(1, 0.5, 0),
  //                            Vec3(-1, -0.5, 0)};
  // std::vector<double> masses = {50, 40};

  // int gridPoints = 64;
  // int dims[] = {gridPoints, gridPoints, gridPoints};
  // KissFFTAdapter<float> fftAdapter(dims, 3);

  // PMMethod pm(gridPoints, fftAdapter);
  // pm.run(state, masses, 15.0, 0.05, 0.5, InterpolationScheme::CIC);
  // ppMethod(state, masses, 10.0);

  glfwInit();
  glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
  glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
  glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);

  GLFWwindow* window = glfwCreateWindow(800, 600, "LearnOpenGL", NULL, NULL);
  if (window == nullptr) {
    std::cerr << "Failed to create GLFW window\n";
    glfwTerminate();
    return -1;
  }
  glfwMakeContextCurrent(window);

  GLenum err = glewInit();
  if (err != GLEW_OK) {
    std::cerr << "Error: " << glewGetErrorString(err) << '\n';
    glfwTerminate();
    return -1;
  }

  glViewport(0, 0, 800, 600);
  glfwSetFramebufferSizeCallback(window, framebufferSizeCallback);

  // render loop
  while (!glfwWindowShouldClose(window)) {
    glClearColor(0.2f, 0.3f, 0.3f, 1.0f);
    glClear(GL_COLOR_BUFFER_BIT);

    glfwSwapBuffers(window);
    glfwPollEvents();
  }

  glfwTerminate();
}

void framebufferSizeCallback(GLFWwindow* window, int width, int height) {
  glViewport(0, 0, width, height);
}