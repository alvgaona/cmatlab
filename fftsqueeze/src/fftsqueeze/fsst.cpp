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
#include "fftsqueeze_data.h"
#include "fsstParser.h"
#include "rt_nonfinite.h"
#include "coder_array.h"
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
  int b_loop_ub;
  int loop_ub;
  int stride_0_0;
  int stride_0_1;
  int stride_1_0;
  int stride_1_1;
  if (in3.size(0) == 1) {
    loop_ub = in2.size(0);
  } else {
    loop_ub = in3.size(0);
  }
  in1.set_size(loop_ub, in1.size(1));
  if (in3.size(1) == 1) {
    b_loop_ub = in2.size(1);
  } else {
    b_loop_ub = in3.size(1);
  }
  in1.set_size(in1.size(0), b_loop_ub);
  stride_0_0 = (in2.size(0) != 1);
  stride_0_1 = (in2.size(1) != 1);
  stride_1_0 = (in3.size(0) != 1);
  stride_1_1 = (in3.size(1) != 1);
  aux_0_1 = 0;
  aux_1_1 = 0;
  for (int i{0}; i < b_loop_ub; i++) {
    for (int i1{0}; i1 < loop_ub; i1++) {
      double ai;
      double ar;
      double bi;
      double d;
      double in2_im;
      int ar_tmp;
      ar_tmp = i1 * stride_0_0;
      ar = in2[ar_tmp + in2.size(0) * aux_0_1].re;
      ai = in2[ar_tmp + in2.size(0) * aux_0_1].im;
      ar_tmp = i1 * stride_1_0;
      d = in3[ar_tmp + in3.size(0) * aux_1_1].re;
      bi = in3[ar_tmp + in3.size(0) * aux_1_1].im;
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
void fsst(const array<double, 2U> &x, double varargin_1,
          const array<double, 1U> &varargin_2, array<creal_T, 2U> &sst,
          array<double, 1U> &f, array<double, 2U> &t)
{
  __m128d r;
  array<creal_T, 2U> Xx;
  array<creal_T, 2U> stftc;
  array<creal_T, 2U> z;
  array<creal_T, 1U> Sxx;
  array<double, 2U> fcorr;
  array<double, 2U> inds;
  array<double, 2U> tout;
  array<double, 2U> xin;
  array<double, 1U> colIdx;
  array<double, 1U> win;
  array<double, 1U> xp;
  double Fs;
  double Nyq;
  double ai;
  double bi;
  double half_res;
  double hopSize;
  double nCol;
  double y_im_tmp;
  int acoef;
  int b_acoef;
  int b_bcoef;
  int bcoef;
  int i;
  int i1;
  int i2;
  int i3;
  int i4;
  int vectorUB;
  Fs = b_signal::internal::fsst::fsstParser(x, varargin_1, varargin_2, win);
  if (x.size(1) - 1 < 0) {
    tout.set_size(tout.size(0), 0);
  } else {
    tout.set_size(1, x.size(1));
    i = x.size(1) - 1;
    for (i1 = 0; i1 <= i; i1++) {
      tout[i1] = i1;
    }
  }
  tout.set_size(1, tout.size(1));
  i = tout.size(1) - 1;
  bcoef = (tout.size(1) / 2) << 1;
  vectorUB = bcoef - 2;
  for (i1 = 0; i1 <= vectorUB; i1 += 2) {
    r = _mm_loadu_pd(&tout[i1]);
    _mm_storeu_pd(&tout[i1], _mm_div_pd(r, _mm_set1_pd(Fs)));
  }
  for (i1 = bcoef; i1 <= i; i1++) {
    tout[i1] = tout[i1] / Fs;
  }
  i1 = static_cast<int>(std::fmod(static_cast<double>(win.size(0)), 2.0));
  if (win.size(0) == 0) {
    i = 0;
  } else {
    i = i1;
  }
  if (i == 1) {
    i = static_cast<int>((static_cast<double>(win.size(0)) - 1.0) / 2.0);
    xp.set_size((i + x.size(1)) + i);
    for (i2 = 0; i2 < i; i2++) {
      xp[i2] = 0.0;
    }
    acoef = x.size(1);
    for (i2 = 0; i2 < acoef; i2++) {
      xp[i2 + i] = x[i2];
    }
    for (i2 = 0; i2 < i; i2++) {
      xp[(i2 + i) + x.size(1)] = 0.0;
    }
  } else {
    i = static_cast<int>(static_cast<double>(win.size(0)) / 2.0);
    acoef = static_cast<int>((static_cast<double>(win.size(0)) - 2.0) / 2.0);
    xp.set_size((i + x.size(1)) + acoef);
    for (i2 = 0; i2 < i; i2++) {
      xp[i2] = 0.0;
    }
    bcoef = x.size(1);
    for (i2 = 0; i2 < bcoef; i2++) {
      xp[i2 + i] = x[i2];
    }
    for (i2 = 0; i2 < acoef; i2++) {
      xp[(i2 + i) + x.size(1)] = 0.0;
    }
  }
  hopSize = static_cast<double>(win.size(0)) -
            (static_cast<double>(win.size(0)) - 1.0);
  nCol = std::trunc((static_cast<double>(xp.size(0)) -
                     (static_cast<double>(win.size(0)) - 1.0)) /
                    hopSize);
  i = win.size(0);
  i2 = static_cast<int>(nCol);
  xin.set_size(win.size(0), i2);
  acoef = win.size(0) * static_cast<int>(nCol);
  for (i3 = 0; i3 < acoef; i3++) {
    xin[i3] = 0.0;
  }
  for (acoef = 0; acoef < i2; acoef++) {
    half_res = hopSize * ((static_cast<double>(acoef) + 1.0) - 1.0);
    if (half_res + 1.0 > static_cast<double>(win.size(0)) + half_res) {
      i3 = 1;
    } else {
      i3 = static_cast<int>(half_res + 1.0);
    }
    for (i4 = 0; i4 < i; i4++) {
      xin[i4 + xin.size(0) * acoef] = xp[(i3 + i4) - 1];
    }
  }
  bsxfun(win, xin, fcorr);
  computeDFTviaFFT(fcorr, static_cast<double>(fcorr.size(0)),
                   static_cast<double>(win.size(0)), Fs, Xx, f);
  b_signal::internal::spectral::dtwin(win, Fs, xp);
  bsxfun(xp, xin, fcorr);
  computeDFTviaFFT(fcorr, static_cast<double>(fcorr.size(0)),
                   static_cast<double>(win.size(0)), Fs, stftc, xp);
  if ((stftc.size(0) == Xx.size(0)) && (stftc.size(1) == Xx.size(1))) {
    fcorr.set_size(stftc.size(0), stftc.size(1));
    acoef = stftc.size(0) * stftc.size(1);
    for (i2 = 0; i2 < acoef; i2++) {
      Nyq = stftc[i2].re;
      ai = stftc[i2].im;
      half_res = Xx[i2].re;
      bi = Xx[i2].im;
      if (bi == 0.0) {
        if (ai == 0.0) {
          nCol = 0.0;
        } else {
          nCol = ai / half_res;
        }
      } else if (half_res == 0.0) {
        if (Nyq == 0.0) {
          nCol = 0.0;
        } else {
          nCol = -(Nyq / bi);
        }
      } else {
        nCol = std::abs(half_res);
        hopSize = std::abs(bi);
        if (nCol > hopSize) {
          nCol = bi / half_res;
          nCol = (ai - nCol * Nyq) / (half_res + nCol * bi);
        } else if (hopSize == nCol) {
          if (half_res > 0.0) {
            half_res = 0.5;
          } else {
            half_res = -0.5;
          }
          if (bi > 0.0) {
            hopSize = 0.5;
          } else {
            hopSize = -0.5;
          }
          nCol = (ai * half_res - Nyq * hopSize) / nCol;
        } else {
          nCol = half_res / bi;
          nCol = (nCol * ai - Nyq) / (bi + nCol * half_res);
        }
      }
      fcorr[i2] = -nCol;
    }
  } else {
    binary_expand_op(fcorr, stftc, Xx);
  }
  acoef = fcorr.size(0) * fcorr.size(1);
  for (vectorUB = 0; vectorUB < acoef; vectorUB++) {
    if (std::isinf(fcorr[vectorUB]) || std::isnan(fcorr[vectorUB])) {
      fcorr[vectorUB] = 0.0;
    }
  }
  i2 = fcorr.size(1);
  xin.set_size(fcorr.size(0), fcorr.size(1));
  for (i3 = 0; i3 < acoef; i3++) {
    xin[i3] = fcorr[i3];
  }
  acoef = fcorr.size(0);
  bcoef = f.size(0);
  if (acoef <= bcoef) {
    bcoef = acoef;
  }
  if (fcorr.size(0) == 1) {
    acoef = f.size(0);
  } else if (f.size(0) == 1) {
    acoef = fcorr.size(0);
  } else if (f.size(0) == fcorr.size(0)) {
    acoef = f.size(0);
  } else {
    acoef = bcoef;
  }
  fcorr.set_size(acoef, i2);
  if ((acoef != 0) && (i2 != 0)) {
    b_bcoef = (xin.size(1) != 1);
    b_acoef = (f.size(0) != 1);
    bcoef = (xin.size(0) != 1);
    for (int k{0}; k < i2; k++) {
      vectorUB = b_bcoef * k;
      for (i = 0; i < acoef; i++) {
        fcorr[i + fcorr.size(0) * k] =
            f[b_acoef * i] + xin[bcoef * i + xin.size(0) * vectorUB];
      }
    }
  }
  acoef = fcorr.size(1);
  bcoef = tout.size(1);
  if (acoef <= bcoef) {
    bcoef = acoef;
  }
  if (fcorr.size(1) == 1) {
    acoef = tout.size(1);
  } else if (tout.size(1) == 1) {
    acoef = fcorr.size(1);
  } else if (tout.size(1) == fcorr.size(1)) {
    acoef = tout.size(1);
  } else {
    acoef = bcoef;
  }
  xin.set_size(fcorr.size(0), acoef);
  if ((fcorr.size(0) != 0) && (acoef != 0)) {
    b_acoef = (tout.size(1) != 1);
    for (int k{0}; k < acoef; k++) {
      bcoef = b_acoef * k;
      i2 = xin.size(0);
      for (i = 0; i < i2; i++) {
        xin[i + xin.size(0) * k] = tout[bcoef];
      }
    }
  }
  if (win.size(0) - 1 < 0) {
    inds.set_size(1, 0);
  } else {
    inds.set_size(1, win.size(0));
    i = win.size(0) - 1;
    for (i2 = 0; i2 <= i; i2++) {
      inds[i2] = i2;
    }
  }
  y_im_tmp = static_cast<double>(win.size(0)) / 2.0;
  nCol = std::floor(y_im_tmp) * -6.2831853071795862;
  i = inds.size(1);
  z.set_size(1, inds.size(1));
  acoef = win.size(0);
  for (i2 = 0; i2 < i; i2++) {
    ai = nCol * inds[i2];
    if (ai == 0.0) {
      z[i2].re = -0.0 / static_cast<double>(acoef);
      z[i2].im = 0.0;
    } else {
      z[i2].re = 0.0;
      z[i2].im = ai / static_cast<double>(acoef);
    }
  }
  for (int k{0}; k < i; k++) {
    if (z[k].re == 0.0) {
      half_res = z[k].im;
      z[k].re = std::cos(half_res);
      z[k].im = std::sin(half_res);
    } else if (z[k].im == 0.0) {
      z[k].re = std::exp(z[k].re);
      z[k].im = 0.0;
    } else if (std::isinf(z[k].im) && std::isinf(z[k].re) && (z[k].re < 0.0)) {
      z[k].re = 0.0;
      z[k].im = 0.0;
    } else {
      nCol = std::exp(z[k].re / 2.0);
      half_res = z[k].im;
      z[k].re = nCol * (nCol * std::cos(half_res));
      z[k].im = nCol * (nCol * std::sin(half_res));
    }
  }
  Sxx.set_size(inds.size(1));
  for (i2 = 0; i2 < i; i2++) {
    Sxx[i2].re = z[i2].re;
    Sxx[i2].im = -z[i2].im;
  }
  acoef = Sxx.size(0);
  bcoef = Xx.size(0);
  if (acoef <= bcoef) {
    bcoef = acoef;
  }
  if (Sxx.size(0) == 1) {
    i2 = Xx.size(0);
  } else if (Xx.size(0) == 1) {
    i2 = Sxx.size(0);
  } else if (Xx.size(0) == Sxx.size(0)) {
    i2 = Xx.size(0);
  } else {
    i2 = bcoef;
  }
  i3 = Xx.size(1);
  stftc.set_size(i2, Xx.size(1));
  if ((i2 != 0) && (Xx.size(1) != 0)) {
    b_acoef = (Xx.size(1) != 1);
    acoef = (Xx.size(0) != 1);
    b_bcoef = (Sxx.size(0) != 1);
    for (int k{0}; k < i3; k++) {
      bcoef = b_acoef * k;
      for (i = 0; i < i2; i++) {
        i4 = acoef * i;
        vectorUB = b_bcoef * i;
        half_res = Xx[i4 + Xx.size(0) * bcoef].re;
        nCol = Sxx[vectorUB].im;
        hopSize = Xx[i4 + Xx.size(0) * bcoef].im;
        Nyq = Sxx[vectorUB].re;
        stftc[i + stftc.size(0) * k].re = half_res * Nyq - hopSize * nCol;
        stftc[i + stftc.size(0) * k].im = half_res * nCol + hopSize * Nyq;
      }
    }
  }
  nCol = f[0];
  hopSize = static_cast<double>(f.size(0)) - 1.0;
  Nyq = f[f.size(0) - 1] - f[0];
  i = fcorr.size(0) * fcorr.size(1);
  xp.set_size(i);
  bcoef = (i / 2) << 1;
  vectorUB = bcoef - 2;
  for (i2 = 0; i2 <= vectorUB; i2 += 2) {
    r = _mm_loadu_pd(&fcorr[i2]);
    _mm_storeu_pd(&xp[i2],
                  _mm_div_pd(_mm_mul_pd(_mm_sub_pd(r, _mm_set1_pd(nCol)),
                                        _mm_set1_pd(hopSize)),
                             _mm_set1_pd(Nyq)));
  }
  for (i2 = bcoef; i2 < i; i2++) {
    xp[i2] = (fcorr[i2] - nCol) * hopSize / Nyq;
  }
  b_acoef = xp.size(0);
  for (int k{0}; k < b_acoef; k++) {
    xp[k] = std::round(xp[k]);
  }
  for (i2 = 0; i2 < b_acoef; i2++) {
    nCol = xp[i2];
    bcoef = f.size(0);
    if (bcoef == 0) {
      if (nCol == 0.0) {
        nCol = 0.0;
      }
    } else if (std::isnan(nCol) || std::isinf(nCol)) {
      nCol = rtNaN;
    } else if (nCol == 0.0) {
      nCol = 0.0;
    } else {
      nCol = std::fmod(nCol, static_cast<double>(bcoef));
      if (nCol == 0.0) {
        nCol = 0.0;
      } else if (nCol < 0.0) {
        nCol += static_cast<double>(bcoef);
      }
    }
    xp[i2] = nCol + 1.0;
  }
  nCol = tout[0];
  hopSize = static_cast<double>(tout.size(1)) - 1.0;
  Nyq = tout[tout.size(1) - 1] - tout[0];
  i = xin.size(0) * xin.size(1);
  colIdx.set_size(i);
  bcoef = (i / 2) << 1;
  vectorUB = bcoef - 2;
  for (i2 = 0; i2 <= vectorUB; i2 += 2) {
    r = _mm_loadu_pd(&xin[i2]);
    _mm_storeu_pd(&colIdx[i2],
                  _mm_div_pd(_mm_mul_pd(_mm_sub_pd(r, _mm_set1_pd(nCol)),
                                        _mm_set1_pd(hopSize)),
                             _mm_set1_pd(Nyq)));
  }
  for (i2 = bcoef; i2 < i; i2++) {
    colIdx[i2] = (xin[i2] - nCol) * hopSize / Nyq;
  }
  acoef = colIdx.size(0);
  for (int k{0}; k < acoef; k++) {
    colIdx[k] = std::round(colIdx[k]);
  }
  bcoef = (colIdx.size(0) / 2) << 1;
  vectorUB = bcoef - 2;
  for (i2 = 0; i2 <= vectorUB; i2 += 2) {
    r = _mm_loadu_pd(&colIdx[i2]);
    _mm_storeu_pd(&colIdx[i2], _mm_add_pd(r, _mm_set1_pd(1.0)));
  }
  for (i2 = bcoef; i2 < acoef; i2++) {
    colIdx[i2] = colIdx[i2] + 1.0;
  }
  i = stftc.size(0) * stftc.size(1);
  Sxx.set_size(i);
  for (i2 = 0; i2 < i; i2++) {
    Sxx[i2] = stftc[i2];
  }
  b_bcoef = tout.size(1);
  stftc.set_size(f.size(0), tout.size(1));
  acoef = f.size(0) * tout.size(1);
  for (i2 = 0; i2 < acoef; i2++) {
    stftc[i2].re = 0.0;
    stftc[i2].im = 0.0;
  }
  for (vectorUB = 0; vectorUB < b_acoef; vectorUB++) {
    if ((xp[vectorUB] >= 1.0) && (xp[vectorUB] <= f.size(0)) &&
        (colIdx[vectorUB] >= 1.0) && (colIdx[vectorUB] <= tout.size(1))) {
      i2 = static_cast<int>(xp[vectorUB]) - 1;
      i3 = static_cast<int>(colIdx[vectorUB]) - 1;
      stftc[i2 + stftc.size(0) * i3].re =
          stftc[i2 + stftc.size(0) * i3].re + Sxx[vectorUB].re;
      stftc[i2 + stftc.size(0) * i3].im =
          stftc[i2 + stftc.size(0) * i3].im + Sxx[vectorUB].im;
    }
  }
  if (std::isnan(Fs)) {
    nCol = 6.2831853071795862;
  } else {
    nCol = Fs;
  }
  hopSize = nCol / static_cast<double>(win.size(0));
  if (win.size(0) - 1 < 0) {
    inds.set_size(inds.size(0), 0);
  } else {
    inds.set_size(1, win.size(0));
    i = win.size(0) - 1;
    for (i2 = 0; i2 <= i; i2++) {
      inds[i2] = i2;
    }
  }
  inds.set_size(1, inds.size(1));
  i = inds.size(1) - 1;
  bcoef = (inds.size(1) / 2) << 1;
  vectorUB = bcoef - 2;
  for (i2 = 0; i2 <= vectorUB; i2 += 2) {
    r = _mm_loadu_pd(&inds[i2]);
    _mm_storeu_pd(&inds[i2], _mm_mul_pd(_mm_set1_pd(hopSize), r));
  }
  for (i2 = bcoef; i2 <= i; i2++) {
    inds[i2] = hopSize * inds[i2];
  }
  Nyq = nCol / 2.0;
  half_res = hopSize / 2.0;
  if (i1 != 0) {
    bi = (static_cast<double>(win.size(0)) + 1.0) / 2.0;
    inds[static_cast<int>(bi) - 1] = Nyq - half_res;
    inds[static_cast<int>(bi)] = Nyq + half_res;
  } else {
    bi = y_im_tmp + 1.0;
    inds[static_cast<int>(y_im_tmp + 1.0) - 1] = Nyq;
  }
  inds[win.size(0) - 1] = nCol - hopSize;
  acoef = static_cast<int>(bi);
  f.set_size(acoef);
  for (i1 = 0; i1 < acoef; i1++) {
    f[i1] = inds[i1];
  }
  for (i1 = 0; i1 < b_bcoef; i1++) {
    for (i2 = 0; i2 < acoef; i2++) {
      stftc[i2 + static_cast<int>(bi) * i1] = stftc[i2 + stftc.size(0) * i1];
    }
  }
  stftc.set_size(acoef, stftc.size(1));
  sst.set_size(acoef, tout.size(1));
  acoef = static_cast<int>(bi) * tout.size(1);
  for (i1 = 0; i1 < acoef; i1++) {
    sst[i1] = stftc[i1];
  }
  t.set_size(1, tout.size(1));
  for (i1 = 0; i1 < b_bcoef; i1++) {
    t[i1] = tout[i1];
  }
}

} // namespace coder

// End of code generation (fsst.cpp)
