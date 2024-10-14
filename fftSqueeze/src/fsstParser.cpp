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
#include "casyi.h"
#include "cmlri.h"
#include "gammaln.h"
#include "log.h"
#include "rt_nonfinite.h"
#include "coder_array.h"
#include <algorithm>
#include <cmath>

// Function Definitions
namespace coder {
namespace b_signal {
namespace internal {
namespace fsst {
double fsstParser(const array<double, 2U> &x, double varargin_1,
                  const double varargin_2[128], double win_data[],
                  int &win_size)
{
  creal_T hz;
  creal_T tmp;
  creal_T zd;
  double Fs;
  double d;
  int N;
  N = static_cast<int>(std::fmin(256.0, static_cast<double>(x.size(1))));
  if (N > 1) {
    int iseven;
    int mid;
    iseven = 1 - static_cast<int>(static_cast<unsigned int>(N) & 1U);
    mid = (N >> 1) + 1;
    if (mid <= N) {
      zd.im = 0.0;
    }
    for (int k{mid}; k <= N; k++) {
      Fs = static_cast<double>(iseven + ((k - mid) << 1)) /
           (static_cast<double>(N) - 1.0);
      Fs = 10.0 * std::sqrt((1.0 - Fs) * (Fs + 1.0));
      zd.re = Fs;
      if (!std::isnan(Fs)) {
        double az;
        boolean_T guard1;
        if (Fs > 0.0) {
          az = Fs;
        } else {
          az = 0.0;
        }
        guard1 = false;
        if (az <= 2.0) {
          int i;
          int nw;
          nw = 0;
          if ((Fs > 0.0) && (!(Fs < 2.2250738585072014E-305))) {
            hz.re = 0.5 * Fs;
            hz.im = 0.0;
            if (Fs > 4.7170688552396617E-153) {
              Fs = hz.re * hz.re;
              if (!(Fs > 0.0)) {
                Fs = 0.0;
              }
            } else {
              Fs = 0.0;
            }
            b_log(hz);
            d = 1.0;
            gammaln(d);
            hz.re = hz.re * 0.0 - d;
            if (!(hz.re > -700.92179369444591)) {
              nw = 1;
              if (Fs > 0.0) {
                nw = -1;
              }
            }
          }
          if (nw < 0) {
            i = 1;
          } else {
            i = nw;
          }
          if ((1 - i != 0) && (nw < 0)) {
            guard1 = true;
          }
        } else {
          guard1 = true;
        }
        if (guard1) {
          if (az < 21.784271729432426) {
            cmlri(zd, tmp);
          } else {
            casyi(zd, tmp);
          }
        }
      }
    }
  }
  Fs = varargin_1;
  win_size = 128;
  std::copy(&varargin_2[0], &varargin_2[128], &win_data[0]);
  return Fs;
}

} // namespace fsst
} // namespace internal
} // namespace b_signal
} // namespace coder

// End of code generation (fsstParser.cpp)
