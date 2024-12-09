#
# Declare kissfft library
#

add_library(kissfft
  kiss_fft.c
  kfc.c
  kiss_fftnd.c
  kiss_fftndr.c
  kiss_fftr.c)

target_include_directories(kissfft PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>)

