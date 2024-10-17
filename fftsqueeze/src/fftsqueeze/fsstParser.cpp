//
// Trial License - for use to evaluate programs for possible purchase as
// an end-user only.
//
// fsstParser.cpp
//
// Code generation for function 'fsstParser'
//

// Include files
#include "fsstParser.h"
#include "besseli.h"
#include "fftsqueeze_data.h"
#include "rt_nonfinite.h"
#include "coder_array.h"
#include <cmath>

// Function Definitions
namespace coder {
namespace b_signal {
namespace internal {
namespace fsst {
double fsstParser(const array<double, 2U> &x, double varargin_1,
                  const array<double, 1U> &varargin_2, array<double, 1U> &win)
{
  creal_T dc;
  double Fs;
  double a;
  double z_im;
  int iseven;
  int mid;
  int nw;
  nw = static_cast<int>(std::fmin(256.0, static_cast<double>(x.size(1))));
  win.set_size(nw);
  if (nw <= 1) {
    win.set_size(nw);
    for (iseven = 0; iseven < nw; iseven++) {
      win[0] = 1.0;
    }
  } else {
    iseven = 1 - static_cast<int>(static_cast<unsigned int>(nw) & 1U);
    mid = (nw >> 1) + 1;
    for (int k{mid}; k <= nw; k++) {
      Fs = static_cast<double>(iseven + ((k - mid) << 1)) /
           (static_cast<double>(nw) - 1.0);
      dc = besseli(10.0 * std::sqrt((1.0 - Fs) * (Fs + 1.0)));
      if (dc.im == 0.0) {
        Fs = dc.re / 2815.7166284662549;
        z_im = 0.0;
      } else if (dc.re == 0.0) {
        Fs = 0.0;
        z_im = dc.im / 2815.7166284662549;
      } else {
        Fs = dc.re / 2815.7166284662549;
        z_im = dc.im / 2815.7166284662549;
      }
      a = std::abs(Fs);
      Fs = std::abs(z_im);
      if (a < Fs) {
        a /= Fs;
        win[k - 1] = Fs * std::sqrt(a * a + 1.0);
      } else if (a > Fs) {
        Fs /= a;
        win[k - 1] = a * std::sqrt(Fs * Fs + 1.0);
      } else if (std::isnan(Fs)) {
        win[k - 1] = rtNaN;
      } else {
        win[k - 1] = a * 1.4142135623730951;
      }
    }
    for (int k{0}; k <= mid - 2; k++) {
      win[k] = win[(nw - k) - 1];
    }
  }
  if (varargin_2.size(0) != 0) {
    if (varargin_2.size(0) == 1) {
      if (varargin_2[0] == std::floor(varargin_2[0])) {
        nw = static_cast<int>(varargin_2[0]);
      } else {
        nw = static_cast<int>(std::round(varargin_2[0]));
      }
      win.set_size(nw);
      if (nw <= 1) {
        win.set_size(nw);
        for (iseven = 0; iseven < nw; iseven++) {
          win[iseven] = 1.0;
        }
      } else {
        iseven = 1 - static_cast<int>(static_cast<unsigned int>(nw) & 1U);
        mid = (nw >> 1) + 1;
        for (int k{mid}; k <= nw; k++) {
          Fs = static_cast<double>(iseven + ((k - mid) << 1)) /
               (static_cast<double>(nw) - 1.0);
          dc = besseli(10.0 * std::sqrt((1.0 - Fs) * (Fs + 1.0)));
          if (dc.im == 0.0) {
            Fs = dc.re / 2815.7166284662549;
            z_im = 0.0;
          } else if (dc.re == 0.0) {
            Fs = 0.0;
            z_im = dc.im / 2815.7166284662549;
          } else {
            Fs = dc.re / 2815.7166284662549;
            z_im = dc.im / 2815.7166284662549;
          }
          a = std::abs(Fs);
          Fs = std::abs(z_im);
          if (a < Fs) {
            a /= Fs;
            win[k - 1] = Fs * std::sqrt(a * a + 1.0);
          } else if (a > Fs) {
            Fs /= a;
            win[k - 1] = a * std::sqrt(Fs * Fs + 1.0);
          } else if (std::isnan(Fs)) {
            win[k - 1] = rtNaN;
          } else {
            win[k - 1] = a * 1.4142135623730951;
          }
        }
        for (int k{0}; k <= mid - 2; k++) {
          win[k] = win[(nw - k) - 1];
        }
      }
    } else {
      nw = varargin_2.size(0);
      win.set_size(varargin_2.size(0));
      for (iseven = 0; iseven < nw; iseven++) {
        win[iseven] = varargin_2[iseven];
      }
    }
  }
  return varargin_1;
}

} // namespace fsst
} // namespace internal
} // namespace b_signal
} // namespace coder

// End of code generation (fsstParser.cpp)
