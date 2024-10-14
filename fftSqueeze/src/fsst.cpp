//
// Trial License - for use to evaluate programs for possible purchase as
// an end-user only.
//
// fsst.cpp
//
// Code generation for function 'fsst'
//

// Include files
#include "fsst.h"
#include "bsxfun.h"
#include "computeDFT.h"
#include "dtwin.h"
#include "fftSqueeze_data.h"
#include "fsstParser.h"
#include "rt_nonfinite.h"
#include "coder_array.h"
#include <algorithm>
#include <cmath>
#include <emmintrin.h>

// Function Declarations
static void binary_expand_op(coder::array<double, 2U> &in1,
                             const coder::array<creal_T, 2U> &in2,
                             const coder::array<creal_T, 2U> &in3);

// Function Definitions
static void binary_expand_op(coder::array<double, 2U> &in1,
                             const coder::array<creal_T, 2U> &in2,
                             const coder::array<creal_T, 2U> &in3)
{
  int aux_0_1;
  int aux_1_1;
  int loop_ub;
  int stride_0_1;
  int stride_1_1;
  in1.set_size(128, in1.size(1));
  if (in3.size(1) == 1) {
    loop_ub = in2.size(1);
  } else {
    loop_ub = in3.size(1);
  }
  in1.set_size(in1.size(0), loop_ub);
  stride_0_1 = (in2.size(1) != 1);
  stride_1_1 = (in3.size(1) != 1);
  aux_0_1 = 0;
  aux_1_1 = 0;
  for (int i{0}; i < loop_ub; i++) {
    for (int i1{0}; i1 < 128; i1++) {
      double ai;
      double ar;
      double bi;
      double d;
      double in2_im;
      ar = in2[i1 + 128 * aux_0_1].re;
      ai = in2[i1 + 128 * aux_0_1].im;
      d = in3[i1 + 128 * aux_1_1].re;
      bi = in3[i1 + 128 * aux_1_1].im;
      if (bi == 0.0) {
        if (ai == 0.0) {
          in2_im = 0.0;
        } else {
          in2_im = ai / d;
        }
      } else if (d == 0.0) {
        if (ar == 0.0) {
          in2_im = 0.0;
        } else {
          in2_im = -(ar / bi);
        }
      } else {
        double d1;
        in2_im = std::abs(d);
        d1 = std::abs(bi);
        if (in2_im > d1) {
          in2_im = bi / d;
          in2_im = (ai - in2_im * ar) / (d + in2_im * bi);
        } else if (d1 == in2_im) {
          if (d > 0.0) {
            d1 = 0.5;
          } else {
            d1 = -0.5;
          }
          if (bi > 0.0) {
            d = 0.5;
          } else {
            d = -0.5;
          }
          in2_im = (ai * d1 - ar * d) / in2_im;
        } else {
          in2_im = d / bi;
          in2_im = (in2_im * ai - ar) / (bi + in2_im * d);
        }
      }
      in1[i1 + in1.size(0) * i] = -in2_im;
    }
    aux_1_1 += stride_1_1;
    aux_0_1 += stride_0_1;
  }
}

namespace coder {
int fsst(const array<double, 2U> &x, double varargin_1,
         const double varargin_2[128], array<creal_T, 2U> &sst, double f_data[],
         array<double, 2U> &t)
{
  __m128d r;
  array<creal_T, 2U> Xx;
  array<creal_T, 2U> c;
  array<creal_T, 2U> sstout;
  array<creal_T, 2U> stftc;
  array<double, 2U> fcorr;
  array<double, 2U> tcorr;
  array<double, 2U> tout;
  array<double, 2U> xin;
  array<double, 1U> colIdx;
  array<double, 1U> xp;
  creal_T ez_data[128];
  double win_data[256];
  double b_f_data[128];
  double tmp_data[128];
  double w1_data[128];
  double Fs;
  double Fs1;
  double ai;
  double ar;
  double freq_res;
  int acoef;
  int bcoef;
  int f_size;
  int i;
  int i1;
  int loop_ub;
  int nx_tmp;
  int vectorUB;
  Fs = b_signal::internal::fsst::fsstParser(x, varargin_1, varargin_2, win_data,
                                            bcoef);
  if (x.size(1) - 1 < 0) {
    tout.set_size(tout.size(0), 0);
  } else {
    tout.set_size(1, x.size(1));
    loop_ub = x.size(1) - 1;
    for (i = 0; i <= loop_ub; i++) {
      tout[i] = i;
    }
  }
  tout.set_size(1, tout.size(1));
  loop_ub = tout.size(1) - 1;
  acoef = (tout.size(1) / 2) << 1;
  vectorUB = acoef - 2;
  for (i = 0; i <= vectorUB; i += 2) {
    r = _mm_loadu_pd(&tout[i]);
    _mm_storeu_pd(&tout[i], _mm_div_pd(r, _mm_set1_pd(Fs)));
  }
  for (i = acoef; i <= loop_ub; i++) {
    tout[i] = tout[i] / Fs;
  }
  i = x.size(1) + 127;
  xp.set_size(x.size(1) + 127);
  for (i1 = 0; i1 < 64; i1++) {
    xp[i1] = 0.0;
  }
  loop_ub = x.size(1);
  for (i1 = 0; i1 < loop_ub; i1++) {
    xp[i1 + 64] = x[i1];
  }
  for (i1 = 0; i1 < 63; i1++) {
    xp[(i1 + x.size(1)) + 64] = 0.0;
  }
  xin.set_size(128, xp.size(0) - 127);
  loop_ub = (xp.size(0) - 127) << 7;
  for (i1 = 0; i1 < loop_ub; i1++) {
    xin[i1] = 0.0;
  }
  for (loop_ub = 0; loop_ub <= i - 128; loop_ub++) {
    for (i1 = 0; i1 < 128; i1++) {
      xin[i1 + 128 * loop_ub] = xp[loop_ub + i1];
    }
  }
  bsxfun(win_data, bcoef, xin, fcorr);
  f_size = computeDFTviaFFT(fcorr, Fs, Xx, w1_data);
  loop_ub = b_signal::internal::spectral::dtwin(win_data, Fs, tmp_data);
  bsxfun(tmp_data, loop_ub, xin, fcorr);
  computeDFTviaFFT(fcorr, Fs, stftc, b_f_data);
  if (stftc.size(1) == Xx.size(1)) {
    fcorr.set_size(128, stftc.size(1));
    loop_ub = stftc.size(1) << 7;
    for (i = 0; i < loop_ub; i++) {
      double bi;
      double d;
      ar = stftc[i].re;
      ai = stftc[i].im;
      d = Xx[i].re;
      bi = Xx[i].im;
      if (bi == 0.0) {
        if (ai == 0.0) {
          Fs1 = 0.0;
        } else {
          Fs1 = ai / d;
        }
      } else if (d == 0.0) {
        if (ar == 0.0) {
          Fs1 = 0.0;
        } else {
          Fs1 = -(ar / bi);
        }
      } else {
        Fs1 = std::abs(d);
        freq_res = std::abs(bi);
        if (Fs1 > freq_res) {
          Fs1 = bi / d;
          Fs1 = (ai - Fs1 * ar) / (d + Fs1 * bi);
        } else if (freq_res == Fs1) {
          if (d > 0.0) {
            freq_res = 0.5;
          } else {
            freq_res = -0.5;
          }
          if (bi > 0.0) {
            d = 0.5;
          } else {
            d = -0.5;
          }
          Fs1 = (ai * freq_res - ar * d) / Fs1;
        } else {
          Fs1 = d / bi;
          Fs1 = (Fs1 * ai - ar) / (bi + Fs1 * d);
        }
      }
      fcorr[i] = -Fs1;
    }
  } else {
    binary_expand_op(fcorr, stftc, Xx);
  }
  loop_ub = fcorr.size(1) << 7;
  for (vectorUB = 0; vectorUB < loop_ub; vectorUB++) {
    if (std::isinf(fcorr[vectorUB]) || std::isnan(fcorr[vectorUB])) {
      fcorr[vectorUB] = 0.0;
    }
  }
  i = fcorr.size(1);
  xin.set_size(fcorr.size(0), fcorr.size(1));
  loop_ub = fcorr.size(0) * fcorr.size(1);
  for (i1 = 0; i1 < loop_ub; i1++) {
    xin[i1] = fcorr[i1];
  }
  if (f_size == 1) {
    i1 = 128;
  } else if (f_size == 128) {
    i1 = 128;
  } else {
    i1 = f_size;
  }
  fcorr.set_size(i1, i);
  if (i != 0) {
    bcoef = (xin.size(1) != 1);
    acoef = (f_size != 1);
    for (int k{0}; k < i; k++) {
      loop_ub = bcoef * k;
      for (vectorUB = 0; vectorUB < i1; vectorUB++) {
        fcorr[vectorUB + fcorr.size(0) * k] =
            w1_data[acoef * vectorUB] + xin[vectorUB + 128 * loop_ub];
      }
    }
  }
  loop_ub = fcorr.size(1);
  bcoef = tout.size(1);
  if (loop_ub <= bcoef) {
    bcoef = loop_ub;
  }
  if (fcorr.size(1) == 1) {
    bcoef = tout.size(1);
  } else if (tout.size(1) == 1) {
    bcoef = fcorr.size(1);
  } else if (tout.size(1) == fcorr.size(1)) {
    bcoef = tout.size(1);
  }
  tcorr.set_size(fcorr.size(0), bcoef);
  if (bcoef != 0) {
    acoef = (tout.size(1) != 1);
    for (int k{0}; k < bcoef; k++) {
      loop_ub = acoef * k;
      i = tcorr.size(0);
      for (vectorUB = 0; vectorUB < i; vectorUB++) {
        tcorr[vectorUB + tcorr.size(0) * k] = tout[loop_ub];
      }
    }
  }
  for (int k{0}; k < 128; k++) {
    ai = -402.12385965949352 * static_cast<double>(k);
    if (ai == 0.0) {
      Fs1 = 0.0;
    } else {
      Fs1 = ai / 128.0;
    }
    ez_data[k].re = std::cos(Fs1);
    ez_data[k].im = -std::sin(Fs1);
  }
  i = Xx.size(1);
  c.set_size(128, Xx.size(1));
  if (Xx.size(1) != 0) {
    acoef = (Xx.size(1) != 1);
    for (int k{0}; k < i; k++) {
      loop_ub = acoef * k;
      for (vectorUB = 0; vectorUB < 128; vectorUB++) {
        Fs1 = Xx[vectorUB + 128 * loop_ub].re;
        freq_res = ez_data[vectorUB].im;
        ar = Xx[vectorUB + 128 * loop_ub].im;
        ai = ez_data[vectorUB].re;
        c[vectorUB + c.size(0) * k].re = Fs1 * ai - ar * freq_res;
        c[vectorUB + c.size(0) * k].im = Fs1 * freq_res + ar * ai;
      }
    }
  }
  Fs1 = w1_data[0];
  freq_res = w1_data[f_size - 1] - w1_data[0];
  loop_ub = fcorr.size(0) * fcorr.size(1);
  xp.set_size(loop_ub);
  acoef = (loop_ub / 2) << 1;
  vectorUB = acoef - 2;
  for (i = 0; i <= vectorUB; i += 2) {
    r = _mm_loadu_pd(&fcorr[i]);
    _mm_storeu_pd(
        &xp[i],
        _mm_div_pd(_mm_mul_pd(_mm_sub_pd(r, _mm_set1_pd(Fs1)),
                              _mm_set1_pd(static_cast<double>(f_size) - 1.0)),
                   _mm_set1_pd(freq_res)));
  }
  for (i = acoef; i < loop_ub; i++) {
    xp[i] = (fcorr[i] - Fs1) * (static_cast<double>(f_size) - 1.0) / freq_res;
  }
  nx_tmp = xp.size(0);
  for (int k{0}; k < nx_tmp; k++) {
    xp[k] = std::round(xp[k]);
  }
  for (i = 0; i < nx_tmp; i++) {
    Fs1 = xp[i];
    if (std::isnan(Fs1) || std::isinf(Fs1)) {
      Fs1 = rtNaN;
    } else if (Fs1 == 0.0) {
      Fs1 = 0.0;
    } else {
      Fs1 = std::fmod(Fs1, static_cast<double>(f_size));
      if (Fs1 == 0.0) {
        Fs1 = 0.0;
      } else if (Fs1 < 0.0) {
        Fs1 += static_cast<double>(f_size);
      }
    }
    xp[i] = Fs1 + 1.0;
  }
  Fs1 = tout[0];
  freq_res = static_cast<double>(tout.size(1)) - 1.0;
  ar = tout[tout.size(1) - 1] - tout[0];
  loop_ub = tcorr.size(0) * tcorr.size(1);
  colIdx.set_size(loop_ub);
  acoef = (loop_ub / 2) << 1;
  vectorUB = acoef - 2;
  for (i = 0; i <= vectorUB; i += 2) {
    r = _mm_loadu_pd(&tcorr[i]);
    _mm_storeu_pd(&colIdx[i],
                  _mm_div_pd(_mm_mul_pd(_mm_sub_pd(r, _mm_set1_pd(Fs1)),
                                        _mm_set1_pd(freq_res)),
                             _mm_set1_pd(ar)));
  }
  for (i = acoef; i < loop_ub; i++) {
    colIdx[i] = (tcorr[i] - Fs1) * freq_res / ar;
  }
  bcoef = colIdx.size(0);
  for (int k{0}; k < bcoef; k++) {
    colIdx[k] = std::round(colIdx[k]);
  }
  acoef = (colIdx.size(0) / 2) << 1;
  vectorUB = acoef - 2;
  for (i = 0; i <= vectorUB; i += 2) {
    r = _mm_loadu_pd(&colIdx[i]);
    _mm_storeu_pd(&colIdx[i], _mm_add_pd(r, _mm_set1_pd(1.0)));
  }
  for (i = acoef; i < bcoef; i++) {
    colIdx[i] = colIdx[i] + 1.0;
  }
  bcoef = tout.size(1);
  sstout.set_size(f_size, tout.size(1));
  loop_ub = f_size * tout.size(1);
  for (i = 0; i < loop_ub; i++) {
    sstout[i].re = 0.0;
    sstout[i].im = 0.0;
  }
  for (vectorUB = 0; vectorUB < nx_tmp; vectorUB++) {
    if ((xp[vectorUB] >= 1.0) && (xp[vectorUB] <= f_size) &&
        (colIdx[vectorUB] >= 1.0) && (colIdx[vectorUB] <= tout.size(1))) {
      i = static_cast<int>(xp[vectorUB]) - 1;
      i1 = static_cast<int>(colIdx[vectorUB]) - 1;
      sstout[i + sstout.size(0) * i1].re =
          sstout[i + sstout.size(0) * i1].re + c[vectorUB].re;
      sstout[i + sstout.size(0) * i1].im =
          sstout[i + sstout.size(0) * i1].im + c[vectorUB].im;
    }
  }
  for (i = 0; i < bcoef; i++) {
    for (i1 = 0; i1 < 65; i1++) {
      sstout[i1 + 65 * i] = sstout[i1 + sstout.size(0) * i];
    }
  }
  sstout.set_size(65, sstout.size(1));
  sst.set_size(65, tout.size(1));
  loop_ub = 65 * tout.size(1);
  for (i = 0; i < loop_ub; i++) {
    sst[i] = sstout[i];
  }
  if (std::isnan(Fs)) {
    Fs1 = 6.2831853071795862;
  } else {
    Fs1 = Fs;
  }
  freq_res = Fs1 / 128.0;
  for (i = 0; i < 128; i++) {
    w1_data[i] = freq_res * static_cast<double>(i);
  }
  w1_data[64] = Fs1 / 2.0;
  w1_data[127] = Fs1 - freq_res;
  f_size = 65;
  std::copy(&w1_data[0], &w1_data[65], &f_data[0]);
  t.set_size(1, tout.size(1));
  for (i = 0; i < bcoef; i++) {
    t[i] = tout[i];
  }
  return f_size;
}

} // namespace coder

// End of code generation (fsst.cpp)
