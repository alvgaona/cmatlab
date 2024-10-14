//
// Trial License - for use to evaluate programs for possible purchase as
// an end-user only.
//
// fftSqueeze.cpp
//
// Code generation for function 'fftSqueeze'
//

// Include files
#include "fftSqueeze.h"
#include "fftSqueeze_data.h"
#include "fftSqueeze_initialize.h"
#include "fsst.h"
#include "rt_nonfinite.h"
#include "coder_array.h"

// Function Definitions
void fftSqueeze(const coder::array<double, 2U> &x, double fs,
                const double window[128], coder::array<creal_T, 2U> &s,
                double f_data[], int f_size[1], coder::array<double, 2U> &t)
{
  if (!isInitialized_fftSqueeze) {
    fftSqueeze_initialize();
  }
  f_size[0] = coder::fsst(x, fs, window, s, f_data, t);
}

// End of code generation (fftSqueeze.cpp)
