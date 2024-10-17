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
#include <cmath>
#include <emmintrin.h>

// Function Declarations
static int div_s32(int numerator, int denominator);

// Function Definitions
static int div_s32(int numerator, int denominator)
{
  int quotient;
  if (denominator == 0) {
    if (numerator >= 0) {
      quotient = MAX_int32_T;
    } else {
      quotient = MIN_int32_T;
    }
  } else {
    unsigned int u;
    unsigned int u1;
    if (numerator < 0) {
      u = ~static_cast<unsigned int>(numerator) + 1U;
    } else {
      u = static_cast<unsigned int>(numerator);
    }
    if (denominator < 0) {
      u1 = ~static_cast<unsigned int>(denominator) + 1U;
    } else {
      u1 = static_cast<unsigned int>(denominator);
    }
    u /= u1;
    if ((numerator < 0) != (denominator < 0)) {
      quotient = -static_cast<int>(u);
    } else {
      quotient = static_cast<int>(u);
    }
  }
  return quotient;
}

namespace coder {
void computeDFTviaFFT(const array<double, 2U> &xin, double nx, double nfft,
                      double Fs, array<creal_T, 2U> &Xx, array<double, 1U> &f)
{
  __m128d r;
  array<double, 2U> costab;
  array<double, 2U> costab1q;
  array<double, 2U> sintab;
  array<double, 2U> sintabinv;
  array<double, 2U> wrappedData;
  array<double, 2U> xw;
  double Fs1;
  double Nyq;
  double d;
  double freq_res;
  double half_res;
  int b_remainder;
  int i;
  int k;
  int loop_ub_tmp;
  int nFullPasses;
  int offset;
  int pow2p;
  if (nx > nfft) {
    loop_ub_tmp = static_cast<int>(nfft);
    i = xin.size(1);
    xw.set_size(loop_ub_tmp, xin.size(1));
    offset = static_cast<int>(nfft) * xin.size(1);
    for (pow2p = 0; pow2p < offset; pow2p++) {
      xw[pow2p] = 0.0;
    }
    for (int j{0}; j < i; j++) {
      if (xin.size(0) == 1) {
        wrappedData.set_size(1, loop_ub_tmp);
        for (pow2p = 0; pow2p < loop_ub_tmp; pow2p++) {
          wrappedData[pow2p] = 0.0;
        }
      } else {
        wrappedData.set_size(loop_ub_tmp, 1);
        for (pow2p = 0; pow2p < loop_ub_tmp; pow2p++) {
          wrappedData[pow2p] = 0.0;
        }
      }
      nFullPasses = div_s32(xin.size(0), static_cast<int>(nfft));
      offset = nFullPasses * static_cast<int>(nfft);
      b_remainder = (xin.size(0) - offset) - 1;
      for (k = 0; k <= b_remainder; k++) {
        wrappedData[k] = xin[(offset + k) + xin.size(0) * j];
      }
      pow2p = b_remainder + 2;
      for (k = pow2p; k <= loop_ub_tmp; k++) {
        wrappedData[k - 1] = 0.0;
      }
      for (int b_j{0}; b_j < nFullPasses; b_j++) {
        offset = b_j * static_cast<int>(nfft);
        pow2p = (static_cast<int>(nfft) / 2) << 1;
        b_remainder = pow2p - 2;
        for (k = 0; k <= b_remainder; k += 2) {
          r = _mm_loadu_pd(&wrappedData[k]);
          _mm_storeu_pd(
              &wrappedData[k],
              _mm_add_pd(r,
                         _mm_loadu_pd(&xin[(offset + k) + xin.size(0) * j])));
        }
        for (k = pow2p; k < loop_ub_tmp; k++) {
          wrappedData[k] = wrappedData[k] + xin[(offset + k) + xin.size(0) * j];
        }
      }
      for (pow2p = 0; pow2p < loop_ub_tmp; pow2p++) {
        xw[pow2p + xw.size(0) * j] = wrappedData[pow2p];
      }
    }
  } else {
    xw.set_size(xin.size(0), xin.size(1));
    loop_ub_tmp = xin.size(0) * xin.size(1);
    for (i = 0; i < loop_ub_tmp; i++) {
      xw[i] = xin[i];
    }
  }
  if ((xw.size(0) == 0) || (xw.size(1) == 0) || (static_cast<int>(nfft) == 0)) {
    Xx.set_size(static_cast<int>(nfft), xw.size(1));
    loop_ub_tmp = static_cast<int>(nfft) * xw.size(1);
    for (i = 0; i < loop_ub_tmp; i++) {
      Xx[i].re = 0.0;
      Xx[i].im = 0.0;
    }
  } else {
    boolean_T useRadix2;
    useRadix2 =
        ((static_cast<int>(nfft) > 0) &&
         ((static_cast<int>(nfft) & (static_cast<int>(nfft) - 1)) == 0));
    pow2p = 1;
    if (useRadix2) {
      offset = static_cast<int>(nfft);
    } else {
      if (static_cast<int>(nfft) > 0) {
        nFullPasses = (static_cast<int>(nfft) + static_cast<int>(nfft)) - 1;
        offset = 31;
        if (nFullPasses <= 1) {
          offset = 0;
        } else {
          boolean_T exitg1;
          b_remainder = 0;
          exitg1 = false;
          while ((!exitg1) && (offset - b_remainder > 1)) {
            k = (b_remainder + offset) >> 1;
            pow2p = 1 << k;
            if (pow2p == nFullPasses) {
              offset = k;
              exitg1 = true;
            } else if (pow2p > nFullPasses) {
              offset = k;
            } else {
              b_remainder = k;
            }
          }
        }
        pow2p = 1 << offset;
      }
      offset = pow2p;
    }
    Fs1 = 6.2831853071795862 / static_cast<double>(offset);
    nFullPasses = offset / 2 / 2;
    costab1q.set_size(1, nFullPasses + 1);
    costab1q[0] = 1.0;
    offset = static_cast<int>(static_cast<unsigned int>(nFullPasses) >> 1) - 1;
    for (k = 0; k <= offset; k++) {
      costab1q[k + 1] = std::cos(Fs1 * (static_cast<double>(k) + 1.0));
    }
    i = offset + 2;
    for (k = i; k < nFullPasses; k++) {
      costab1q[k] = std::sin(Fs1 * static_cast<double>(nFullPasses - k));
    }
    costab1q[nFullPasses] = 0.0;
    if (!useRadix2) {
      nFullPasses = costab1q.size(1) - 1;
      offset = (costab1q.size(1) - 1) << 1;
      costab.set_size(1, offset + 1);
      sintab.set_size(1, offset + 1);
      costab[0] = 1.0;
      sintab[0] = 0.0;
      sintabinv.set_size(1, offset + 1);
      for (k = 0; k < nFullPasses; k++) {
        sintabinv[k + 1] = costab1q[(nFullPasses - k) - 1];
      }
      i = costab1q.size(1);
      for (k = i; k <= offset; k++) {
        sintabinv[k] = costab1q[k - nFullPasses];
      }
      for (k = 0; k < nFullPasses; k++) {
        costab[k + 1] = costab1q[k + 1];
        sintab[k + 1] = -costab1q[(nFullPasses - k) - 1];
      }
      for (k = i; k <= offset; k++) {
        costab[k] = -costab1q[offset - k];
        sintab[k] = -costab1q[k - nFullPasses];
      }
    } else {
      nFullPasses = costab1q.size(1) - 1;
      offset = (costab1q.size(1) - 1) << 1;
      costab.set_size(1, offset + 1);
      sintab.set_size(1, offset + 1);
      costab[0] = 1.0;
      sintab[0] = 0.0;
      for (k = 0; k < nFullPasses; k++) {
        costab[k + 1] = costab1q[k + 1];
        sintab[k + 1] = -costab1q[(nFullPasses - k) - 1];
      }
      i = costab1q.size(1);
      for (k = i; k <= offset; k++) {
        costab[k] = -costab1q[offset - k];
        sintab[k] = -costab1q[k - nFullPasses];
      }
      sintabinv.set_size(1, 0);
    }
    if (useRadix2) {
      internal::fft::FFTImplementationCallback::r2br_r2dit_trig(
          xw, static_cast<int>(nfft), costab, sintab, Xx);
    } else {
      internal::fft::FFTImplementationCallback::dobluesteinfft(
          xw, pow2p, static_cast<int>(nfft), costab, sintab, sintabinv, Xx);
    }
  }
  if (std::isnan(Fs)) {
    Fs1 = 6.2831853071795862;
  } else {
    Fs1 = Fs;
  }
  freq_res = Fs1 / nfft;
  if (std::isnan(nfft - 1.0)) {
    costab1q.set_size(1, 1);
    costab1q[0] = rtNaN;
  } else if (nfft - 1.0 < 0.0) {
    costab1q.set_size(costab1q.size(0), 0);
  } else {
    costab1q.set_size(1, static_cast<int>(nfft - 1.0) + 1);
    offset = static_cast<int>(nfft - 1.0);
    for (i = 0; i <= offset; i++) {
      costab1q[i] = i;
    }
  }
  costab1q.set_size(1, costab1q.size(1));
  offset = costab1q.size(1) - 1;
  pow2p = (costab1q.size(1) / 2) << 1;
  b_remainder = pow2p - 2;
  for (i = 0; i <= b_remainder; i += 2) {
    r = _mm_loadu_pd(&costab1q[i]);
    _mm_storeu_pd(&costab1q[i], _mm_mul_pd(_mm_set1_pd(freq_res), r));
  }
  for (i = pow2p; i <= offset; i++) {
    costab1q[i] = freq_res * costab1q[i];
  }
  Nyq = Fs1 / 2.0;
  half_res = freq_res / 2.0;
  if (std::isnan(nfft) || std::isinf(nfft)) {
    d = rtNaN;
  } else {
    d = std::fmod(nfft, 2.0);
  }
  if (d != 0.0) {
    d = (nfft + 1.0) / 2.0;
    costab1q[static_cast<int>(d) - 1] = Nyq - half_res;
    costab1q[static_cast<int>(static_cast<unsigned int>(d))] = Nyq + half_res;
  } else {
    costab1q[static_cast<int>(nfft / 2.0 + 1.0) - 1] = Nyq;
  }
  costab1q[static_cast<int>(nfft) - 1] = Fs1 - freq_res;
  offset = costab1q.size(1);
  f.set_size(costab1q.size(1));
  for (i = 0; i < offset; i++) {
    f[i] = costab1q[i];
  }
}

} // namespace coder

// End of code generation (computeDFT.cpp)
