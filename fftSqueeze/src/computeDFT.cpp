//
// Trial License - for use to evaluate programs for possible purchase as
// an end-user only.
//
// computeDFT.cpp
//
// Code generation for function 'computeDFT'
//

// Include files
#include "computeDFT.h"
#include "FFTImplementationCallback.h"
#include "rt_nonfinite.h"
#include "coder_array.h"
#include <algorithm>
#include <cmath>

// Function Definitions
namespace coder {
int computeDFTviaFFT(const array<double, 2U> &xin, double Fs,
                     array<creal_T, 2U> &Xx, double f_data[])
{
  double w1_data[128];
  double Fs1;
  double freq_res;
  int f_size;
  if (xin.size(1) == 0) {
    Xx.set_size(128, 0);
  } else {
    double costab_data[65];
    double sintab_data[65];
    double costab1q_data[33];
    costab1q_data[0] = 1.0;
    for (int k{0}; k < 16; k++) {
      costab1q_data[k + 1] =
          std::cos(0.049087385212340517 * (static_cast<double>(k) + 1.0));
    }
    for (int k{0}; k < 15; k++) {
      costab1q_data[k + 17] = std::sin(
          0.049087385212340517 * (32.0 - (static_cast<double>(k) + 17.0)));
    }
    costab1q_data[32] = 0.0;
    costab_data[0] = 1.0;
    sintab_data[0] = 0.0;
    for (int k{0}; k < 32; k++) {
      Fs1 = costab1q_data[k + 1];
      costab_data[k + 1] = Fs1;
      freq_res = -costab1q_data[31 - k];
      sintab_data[k + 1] = freq_res;
      costab_data[k + 33] = freq_res;
      sintab_data[k + 33] = -Fs1;
    }
    internal::fft::FFTImplementationCallback::r2br_r2dit_trig(xin, costab_data,
                                                              sintab_data, Xx);
  }
  if (std::isnan(Fs)) {
    Fs1 = 6.2831853071795862;
  } else {
    Fs1 = Fs;
  }
  freq_res = Fs1 / 128.0;
  for (int k{0}; k < 128; k++) {
    w1_data[k] = freq_res * static_cast<double>(k);
  }
  w1_data[64] = Fs1 / 2.0;
  w1_data[127] = Fs1 - freq_res;
  f_size = 128;
  std::copy(&w1_data[0], &w1_data[128], &f_data[0]);
  return f_size;
}

} // namespace coder

// End of code generation (computeDFT.cpp)
