#define cudaCheckErrors(msg)                                                                       \
  do {                                                                                             \
    cudaError_t __err = cudaGetLastError();                                                        \
    if (__err != cudaSuccess) {                                                                    \
      fprintf(stderr, "Fatal error: %s (%s at %s:%d)\n", msg, cudaGetErrorString(__err), __FILE__, \
              __LINE__);                                                                           \
      fprintf(stderr, "*** FAILED - ABORTING\n");                                                  \
      exit(1);                                                                                     \
    }                                                                                              \
  } while (0)

#define cudaTime(event, ...)       \
  do {                             \
    cudaEventRecord(event##Start); \
    __VA_ARGS__;                   \
    cudaEventRecord(event##Stop);  \
  } while (0)

#define cudaAccTime(acc, event)                            \
  do {                                                     \
    acc = 0;                                               \
    cudaEventElapsedTime(&acc, event##Start, event##Stop); \
    event##Ms += acc;                                      \
  } while (0)

#define hostTime(event, ...)                               \
  do {                                                     \
    auto event##Start_ = std::chrono::steady_clock::now(); \
    __VA_ARGS__;                                           \
    auto event##Stop_ = std::chrono::steady_clock::now();  \
    event##Time##Ms += toMs(event##Stop_ - event##Start_); \
  } while (0)