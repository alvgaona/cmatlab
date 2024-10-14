//
// Trial License - for use to evaluate programs for possible purchase as
// an end-user only.
//
// FFTImplementationCallback.cpp
//
// Code generation for function 'FFTImplementationCallback'
//

// Include files
#include "FFTImplementationCallback.h"
#include "rt_nonfinite.h"
#include "coder_array.h"
#include "omp.h"
#include <cstring>

// Function Definitions
namespace coder {
namespace internal {
namespace fft {
void FFTImplementationCallback::doHalfLengthRadix2(const array<double, 2U> &x,
                                                   int xoffInit,
                                                   creal_T y_data[],
                                                   const double costab_data[],
                                                   const double sintab_data[])
{
  creal_T reconVar1_data[64];
  creal_T reconVar2_data[64];
  double hcostab_data[32];
  double hsintab_data[32];
  double b_temp_re_tmp;
  double im;
  double re;
  double temp2_im;
  double temp2_re;
  double temp_im;
  double temp_im_tmp;
  double temp_re;
  int bitrevIndex_data[64];
  int wrapIndex_data[64];
  int i;
  int iheight;
  int iy;
  int ju;
  int k;
  boolean_T tst;
  for (i = 0; i < 32; i++) {
    iy = ((i + 1) << 1) - 2;
    hcostab_data[i] = costab_data[iy];
    hsintab_data[i] = sintab_data[iy];
  }
  ju = 0;
  iy = 1;
  for (i = 0; i < 64; i++) {
    temp2_re = sintab_data[i];
    reconVar1_data[i].re = temp2_re + 1.0;
    temp2_im = costab_data[i];
    reconVar1_data[i].im = -temp2_im;
    reconVar2_data[i].re = 1.0 - temp2_re;
    reconVar2_data[i].im = temp2_im;
    if (i + 1 != 1) {
      wrapIndex_data[i] = 65 - i;
    } else {
      wrapIndex_data[0] = 1;
    }
    bitrevIndex_data[i] = 0;
  }
  for (k = 0; k < 63; k++) {
    bitrevIndex_data[k] = iy;
    iy = 64;
    tst = true;
    while (tst) {
      iy >>= 1;
      ju ^= iy;
      tst = ((ju & iy) == 0);
    }
    iy = ju + 1;
  }
  bitrevIndex_data[63] = iy;
  if ((static_cast<unsigned int>(x.size(0)) & 1U) == 0U) {
    tst = true;
    iy = x.size(0);
  } else if (x.size(0) >= 128) {
    tst = true;
    iy = 128;
  } else {
    tst = false;
    iy = x.size(0) - 1;
  }
  ju = static_cast<unsigned char>(iy) >> 1;
  for (i = 0; i < ju; i++) {
    iy = xoffInit + (i << 1);
    k = bitrevIndex_data[i];
    y_data[k - 1].re = x[iy];
    y_data[k - 1].im = x[iy + 1];
  }
  if (!tst) {
    if (ju - 1 < 0) {
      i = xoffInit;
    } else {
      i = xoffInit + (ju << 1);
    }
    y_data[bitrevIndex_data[ju] - 1].re = x[i];
    y_data[bitrevIndex_data[ju] - 1].im = 0.0;
  }
  for (i = 0; i <= 62; i += 2) {
    temp2_re = y_data[i + 1].re;
    temp2_im = y_data[i + 1].im;
    temp_re = temp2_re;
    temp_im = temp2_im;
    re = y_data[i].re;
    im = y_data[i].im;
    temp2_re = re - temp2_re;
    temp2_im = im - temp2_im;
    y_data[i + 1].re = temp2_re;
    y_data[i + 1].im = temp2_im;
    re += temp_re;
    im += temp_im;
    y_data[i].re = re;
    y_data[i].im = im;
  }
  iy = 2;
  ju = 4;
  k = 16;
  iheight = 61;
  while (k > 0) {
    int istart;
    int temp_re_tmp;
    for (i = 0; i < iheight; i += ju) {
      temp_re_tmp = i + iy;
      temp_re = y_data[temp_re_tmp].re;
      temp_im = y_data[temp_re_tmp].im;
      y_data[temp_re_tmp].re = y_data[i].re - temp_re;
      y_data[temp_re_tmp].im = y_data[i].im - temp_im;
      y_data[i].re += temp_re;
      y_data[i].im += temp_im;
    }
    istart = 1;
    for (int j{k}; j < 32; j += k) {
      int ihi;
      temp2_re = hcostab_data[j];
      temp2_im = hsintab_data[j];
      i = istart;
      ihi = istart + iheight;
      while (i < ihi) {
        temp_re_tmp = i + iy;
        b_temp_re_tmp = y_data[temp_re_tmp].im;
        temp_im = y_data[temp_re_tmp].re;
        temp_re = temp2_re * temp_im - temp2_im * b_temp_re_tmp;
        temp_im = temp2_re * b_temp_re_tmp + temp2_im * temp_im;
        y_data[temp_re_tmp].re = y_data[i].re - temp_re;
        y_data[temp_re_tmp].im = y_data[i].im - temp_im;
        y_data[i].re += temp_re;
        y_data[i].im += temp_im;
        i += ju;
      }
      istart++;
    }
    k /= 2;
    iy = ju;
    ju += ju;
    iheight -= iy;
  }
  temp2_re = y_data[0].re;
  temp_im_tmp = y_data[0].im;
  temp_im = temp2_re * reconVar1_data[0].re;
  re = temp2_re * reconVar1_data[0].im;
  temp_re = temp2_re * reconVar2_data[0].re;
  temp2_im = temp2_re * reconVar2_data[0].im;
  y_data[0].re = 0.5 * ((temp_im - temp_im_tmp * reconVar1_data[0].im) +
                        (temp_re - -temp_im_tmp * reconVar2_data[0].im));
  y_data[0].im = 0.5 * ((re + temp_im_tmp * reconVar1_data[0].re) +
                        (temp2_im + -temp_im_tmp * reconVar2_data[0].re));
  y_data[64].re = 0.5 * ((temp_re - temp_im_tmp * reconVar2_data[0].im) +
                         (temp_im - -temp_im_tmp * reconVar1_data[0].im));
  y_data[64].im = 0.5 * ((temp2_im + temp_im_tmp * reconVar2_data[0].re) +
                         (re + -temp_im_tmp * reconVar1_data[0].re));
  for (i = 0; i < 31; i++) {
    double temp2_im_tmp;
    b_temp_re_tmp = y_data[i + 1].re;
    temp_im_tmp = y_data[i + 1].im;
    iy = wrapIndex_data[i + 1];
    temp2_im = y_data[iy - 1].re;
    temp2_im_tmp = y_data[iy - 1].im;
    temp_im = reconVar1_data[i + 1].im;
    temp_re = reconVar1_data[i + 1].re;
    re = reconVar2_data[i + 1].im;
    im = reconVar2_data[i + 1].re;
    y_data[i + 1].re =
        0.5 * ((b_temp_re_tmp * temp_re - temp_im_tmp * temp_im) +
               (temp2_im * im - -temp2_im_tmp * re));
    y_data[i + 1].im =
        0.5 * ((b_temp_re_tmp * temp_im + temp_im_tmp * temp_re) +
               (temp2_im * re + -temp2_im_tmp * im));
    y_data[i + 65].re = 0.5 * ((b_temp_re_tmp * im - temp_im_tmp * re) +
                               (temp2_im * temp_re - -temp2_im_tmp * temp_im));
    y_data[i + 65].im = 0.5 * ((b_temp_re_tmp * re + temp_im_tmp * im) +
                               (temp2_im * temp_im + -temp2_im_tmp * temp_re));
    re = reconVar1_data[iy - 1].im;
    im = reconVar1_data[iy - 1].re;
    temp_im = reconVar2_data[iy - 1].im;
    temp2_re = reconVar2_data[iy - 1].re;
    y_data[iy - 1].re =
        0.5 * ((temp2_im * im - temp2_im_tmp * re) +
               (b_temp_re_tmp * temp2_re - -temp_im_tmp * temp_im));
    y_data[iy - 1].im =
        0.5 * ((temp2_im * re + temp2_im_tmp * im) +
               (b_temp_re_tmp * temp_im + -temp_im_tmp * temp2_re));
    y_data[iy + 63].re = 0.5 * ((temp2_im * temp2_re - temp2_im_tmp * temp_im) +
                                (b_temp_re_tmp * im - -temp_im_tmp * re));
    y_data[iy + 63].im = 0.5 * ((temp2_im * temp_im + temp2_im_tmp * temp2_re) +
                                (b_temp_re_tmp * re + -temp_im_tmp * im));
  }
  temp2_re = y_data[32].re;
  temp_im_tmp = y_data[32].im;
  temp_im = temp2_re * reconVar1_data[32].re;
  re = temp2_re * reconVar1_data[32].im;
  temp_re = temp2_re * reconVar2_data[32].re;
  temp2_im = temp2_re * reconVar2_data[32].im;
  y_data[32].re = 0.5 * ((temp_im - temp_im_tmp * reconVar1_data[32].im) +
                         (temp_re - -temp_im_tmp * reconVar2_data[32].im));
  y_data[32].im = 0.5 * ((re + temp_im_tmp * reconVar1_data[32].re) +
                         (temp2_im + -temp_im_tmp * reconVar2_data[32].re));
  y_data[96].re = 0.5 * ((temp_re - temp_im_tmp * reconVar2_data[32].im) +
                         (temp_im - -temp_im_tmp * reconVar1_data[32].im));
  y_data[96].im = 0.5 * ((temp2_im + temp_im_tmp * reconVar2_data[32].re) +
                         (re + -temp_im_tmp * reconVar1_data[32].re));
}

void FFTImplementationCallback::r2br_r2dit_trig(const array<double, 2U> &x,
                                                const double costab_data[],
                                                const double sintab_data[],
                                                array<creal_T, 2U> &y)
{
  creal_T tmp_data[128];
  int loop_ub;
  int nrows;
  int xoff;
  nrows = x.size(0);
  y.set_size(128, x.size(1));
  if (x.size(0) < 128) {
    y.set_size(128, x.size(1));
    loop_ub = 128 * x.size(1);
    for (int i{0}; i < loop_ub; i++) {
      y[i].re = 0.0;
      y[i].im = 0.0;
    }
  }
  loop_ub = x.size(1);
#pragma omp parallel for num_threads(omp_get_max_threads()) private(           \
        xoff, tmp_data)

  for (int chan = 0; chan < loop_ub; chan++) {
    xoff = chan * nrows;
    if (x.size(0) < 128) {
      std::memset(&tmp_data[0], 0, 128U * sizeof(creal_T));
    }
    FFTImplementationCallback::doHalfLengthRadix2(x, xoff, tmp_data,
                                                  costab_data, sintab_data);
    for (xoff = 0; xoff < 128; xoff++) {
      y[xoff + 128 * chan] = tmp_data[xoff];
    }
  }
}

} // namespace fft
} // namespace internal
} // namespace coder

// End of code generation (FFTImplementationCallback.cpp)
