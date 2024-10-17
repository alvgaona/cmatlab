//
// Trial License - for use to evaluate programs for possible purchase as
// an end-user only.
//
// fftsqueeze.cpp
//
// Code generation for function 'fftsqueeze'
//

// Include files
#include "fftsqueeze.h"
#include "fftsqueeze_data.h"
#include "fftsqueeze_initialize.h"
#include "fsst.h"
#include "rt_nonfinite.h"
#include "coder_array.h"

// Function Definitions
void fftsqueeze(const coder::array<double, 2U> &x, double fs,
                const coder::array<double, 1U> &window,
                coder::array<creal_T, 2U> &s, coder::array<double, 1U> &f,
                coder::array<double, 2U> &t)
{
  if (!isInitialized_fftsqueeze) {
    fftsqueeze_initialize();
  }
  // codegen
  coder::fsst(x, fs, window, s, f, t);
}

// End of code generation (fftsqueeze.cpp)
