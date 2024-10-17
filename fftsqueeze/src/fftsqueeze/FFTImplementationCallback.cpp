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
#include <cmath>

// Function Definitions
namespace coder {
namespace internal {
namespace fft {
void FFTImplementationCallback::doHalfLengthBluestein(
    const array<double, 2U> &x, int xoffInit, array<creal_T, 1U> &y, int nrowsx,
    int nRows, int nfft, const array<creal_T, 1U> &wwc,
    const array<double, 2U> &costab, const array<double, 2U> &sintab,
    const array<double, 2U> &costabinv, const array<double, 2U> &sintabinv)
{
  array<creal_T, 1U> fv;
  array<creal_T, 1U> fy;
  array<creal_T, 1U> reconVar1;
  array<creal_T, 1U> reconVar2;
  array<creal_T, 1U> ytmp;
  array<double, 2U> b_costab;
  array<double, 2U> b_sintab;
  array<double, 2U> costab1q;
  array<double, 2U> hcostabinv;
  array<double, 2U> hsintab;
  array<double, 2U> hsintabinv;
  array<int, 2U> wrapIndex;
  double b_temp_re_tmp;
  double e;
  double temp_im;
  double temp_re;
  double twid_im;
  double twid_re;
  int hnRows;
  int i;
  int istart;
  int j;
  int ju;
  int k;
  int minHnrowsNxBy2;
  int nRowsD2;
  int nd2;
  int nfft_tmp;
  int temp_re_tmp;
  boolean_T tst;
  hnRows = nRows / 2;
  ytmp.set_size(hnRows);
  if (hnRows > nrowsx) {
    ytmp.set_size(hnRows);
    for (minHnrowsNxBy2 = 0; minHnrowsNxBy2 < hnRows; minHnrowsNxBy2++) {
      ytmp[minHnrowsNxBy2].re = 0.0;
      ytmp[minHnrowsNxBy2].im = 0.0;
    }
  }
  if ((static_cast<unsigned int>(x.size(0)) & 1U) == 0U) {
    tst = true;
    j = x.size(0);
  } else if (x.size(0) >= nRows) {
    tst = true;
    j = nRows;
  } else {
    tst = false;
    j = x.size(0) - 1;
  }
  nd2 = nRows << 1;
  e = 6.2831853071795862 / static_cast<double>(nd2);
  istart = nd2 / 2 / 2;
  costab1q.set_size(1, istart + 1);
  costab1q[0] = 1.0;
  nd2 = static_cast<int>(static_cast<unsigned int>(istart) >> 1) - 1;
  for (k = 0; k <= nd2; k++) {
    costab1q[k + 1] = std::cos(e * (static_cast<double>(k) + 1.0));
  }
  minHnrowsNxBy2 = nd2 + 2;
  for (k = minHnrowsNxBy2; k < istart; k++) {
    costab1q[k] = std::sin(e * static_cast<double>(istart - k));
  }
  costab1q[istart] = 0.0;
  istart = costab1q.size(1) - 1;
  nd2 = (costab1q.size(1) - 1) << 1;
  b_costab.set_size(1, nd2 + 1);
  b_sintab.set_size(1, nd2 + 1);
  b_costab[0] = 1.0;
  b_sintab[0] = 0.0;
  for (k = 0; k < istart; k++) {
    b_costab[k + 1] = costab1q[k + 1];
    b_sintab[k + 1] = -costab1q[(istart - k) - 1];
  }
  minHnrowsNxBy2 = costab1q.size(1);
  for (k = minHnrowsNxBy2; k <= nd2; k++) {
    b_costab[k] = -costab1q[nd2 - k];
    b_sintab[k] = -costab1q[k - istart];
  }
  nd2 = static_cast<int>(static_cast<unsigned int>(costab.size(1)) >> 1);
  costab1q.set_size(1, nd2);
  hsintab.set_size(1, nd2);
  hcostabinv.set_size(1, nd2);
  hsintabinv.set_size(1, nd2);
  for (i = 0; i < nd2; i++) {
    minHnrowsNxBy2 = ((i + 1) << 1) - 2;
    costab1q[i] = costab[minHnrowsNxBy2];
    hsintab[i] = sintab[minHnrowsNxBy2];
    hcostabinv[i] = costabinv[minHnrowsNxBy2];
    hsintabinv[i] = sintabinv[minHnrowsNxBy2];
  }
  reconVar1.set_size(hnRows);
  reconVar2.set_size(hnRows);
  wrapIndex.set_size(1, hnRows);
  for (i = 0; i < hnRows; i++) {
    minHnrowsNxBy2 = i << 1;
    e = b_sintab[minHnrowsNxBy2];
    temp_re = b_costab[minHnrowsNxBy2];
    reconVar1[i].re = e + 1.0;
    reconVar1[i].im = -temp_re;
    reconVar2[i].re = 1.0 - e;
    reconVar2[i].im = temp_re;
    if (i + 1 != 1) {
      wrapIndex[i] = (hnRows - i) + 1;
    } else {
      wrapIndex[0] = 1;
    }
  }
  if (j > nRows) {
    j = nRows;
  }
  minHnrowsNxBy2 = j / 2 - 1;
  for (ju = 0; ju <= minHnrowsNxBy2; ju++) {
    temp_re_tmp = (hnRows + ju) - 1;
    temp_re = wwc[temp_re_tmp].re;
    temp_im = wwc[temp_re_tmp].im;
    nd2 = xoffInit + (ju << 1);
    twid_re = x[nd2];
    twid_im = x[nd2 + 1];
    ytmp[ju].re = temp_re * twid_re + temp_im * twid_im;
    ytmp[ju].im = temp_re * twid_im - temp_im * twid_re;
  }
  if (!tst) {
    temp_re_tmp = hnRows + minHnrowsNxBy2;
    temp_re = wwc[temp_re_tmp].re;
    temp_im = wwc[temp_re_tmp].im;
    if (minHnrowsNxBy2 < 0) {
      j = xoffInit;
    } else {
      j = xoffInit + ((minHnrowsNxBy2 + 1) << 1);
    }
    twid_re = x[j];
    ytmp[minHnrowsNxBy2 + 1].re = temp_re * twid_re + temp_im * 0.0;
    ytmp[minHnrowsNxBy2 + 1].im = temp_re * 0.0 - temp_im * twid_re;
    if (minHnrowsNxBy2 + 3 <= hnRows) {
      minHnrowsNxBy2 += 3;
      for (i = minHnrowsNxBy2; i <= hnRows; i++) {
        ytmp[i - 1].re = 0.0;
        ytmp[i - 1].im = 0.0;
      }
    }
  } else if (minHnrowsNxBy2 + 2 <= hnRows) {
    minHnrowsNxBy2 += 2;
    for (i = minHnrowsNxBy2; i <= hnRows; i++) {
      ytmp[i - 1].re = 0.0;
      ytmp[i - 1].im = 0.0;
    }
  }
  nfft_tmp = nfft / 2;
  fy.set_size(nfft_tmp);
  if (nfft_tmp > ytmp.size(0)) {
    fy.set_size(nfft_tmp);
    for (minHnrowsNxBy2 = 0; minHnrowsNxBy2 < nfft_tmp; minHnrowsNxBy2++) {
      fy[minHnrowsNxBy2].re = 0.0;
      fy[minHnrowsNxBy2].im = 0.0;
    }
  }
  j = ytmp.size(0);
  if (j > nfft_tmp) {
    j = nfft_tmp;
  }
  minHnrowsNxBy2 = nfft_tmp - 2;
  nRowsD2 = nfft_tmp / 2;
  k = nRowsD2 / 2;
  nd2 = 0;
  ju = 0;
  for (i = 0; i <= j - 2; i++) {
    fy[nd2] = ytmp[i];
    istart = nfft_tmp;
    tst = true;
    while (tst) {
      istart >>= 1;
      ju ^= istart;
      tst = ((ju & istart) == 0);
    }
    nd2 = ju;
  }
  if (j - 2 < 0) {
    j = 0;
  } else {
    j--;
  }
  fy[nd2] = ytmp[j];
  if (nfft_tmp > 1) {
    for (i = 0; i <= minHnrowsNxBy2; i += 2) {
      b_temp_re_tmp = fy[i + 1].re;
      temp_im = fy[i + 1].im;
      twid_im = fy[i].re;
      e = fy[i].im;
      fy[i + 1].re = twid_im - b_temp_re_tmp;
      fy[i + 1].im = fy[i].im - fy[i + 1].im;
      e += temp_im;
      fy[i].re = twid_im + b_temp_re_tmp;
      fy[i].im = e;
    }
  }
  nd2 = 2;
  minHnrowsNxBy2 = 4;
  ju = ((k - 1) << 2) + 1;
  while (k > 0) {
    for (i = 0; i < ju; i += minHnrowsNxBy2) {
      temp_re_tmp = i + nd2;
      temp_re = fy[temp_re_tmp].re;
      temp_im = fy[temp_re_tmp].im;
      fy[temp_re_tmp].re = fy[i].re - temp_re;
      fy[temp_re_tmp].im = fy[i].im - temp_im;
      fy[i].re = fy[i].re + temp_re;
      fy[i].im = fy[i].im + temp_im;
    }
    istart = 1;
    for (j = k; j < nRowsD2; j += k) {
      int ihi;
      twid_re = costab1q[j];
      twid_im = hsintab[j];
      i = istart;
      ihi = istart + ju;
      while (i < ihi) {
        temp_re_tmp = i + nd2;
        b_temp_re_tmp = fy[temp_re_tmp].im;
        e = fy[temp_re_tmp].re;
        temp_re = twid_re * e - twid_im * b_temp_re_tmp;
        temp_im = twid_re * b_temp_re_tmp + twid_im * e;
        fy[temp_re_tmp].re = fy[i].re - temp_re;
        fy[temp_re_tmp].im = fy[i].im - temp_im;
        fy[i].re = fy[i].re + temp_re;
        fy[i].im = fy[i].im + temp_im;
        i += minHnrowsNxBy2;
      }
      istart++;
    }
    k /= 2;
    nd2 = minHnrowsNxBy2;
    minHnrowsNxBy2 += minHnrowsNxBy2;
    ju -= nd2;
  }
  FFTImplementationCallback::r2br_r2dit_trig_impl(wwc, nfft_tmp, costab1q,
                                                  hsintab, fv);
  for (minHnrowsNxBy2 = 0; minHnrowsNxBy2 < nfft_tmp; minHnrowsNxBy2++) {
    twid_im = fy[minHnrowsNxBy2].re;
    e = fv[minHnrowsNxBy2].im;
    temp_re = fy[minHnrowsNxBy2].im;
    twid_re = fv[minHnrowsNxBy2].re;
    fy[minHnrowsNxBy2].re = twid_im * twid_re - temp_re * e;
    fy[minHnrowsNxBy2].im = twid_im * e + temp_re * twid_re;
  }
  FFTImplementationCallback::r2br_r2dit_trig_impl(fy, nfft_tmp, hcostabinv,
                                                  hsintabinv, fv);
  if (fv.size(0) > 1) {
    e = 1.0 / static_cast<double>(fv.size(0));
    nd2 = fv.size(0);
    for (minHnrowsNxBy2 = 0; minHnrowsNxBy2 < nd2; minHnrowsNxBy2++) {
      fv[minHnrowsNxBy2].re = e * fv[minHnrowsNxBy2].re;
      fv[minHnrowsNxBy2].im = e * fv[minHnrowsNxBy2].im;
    }
  }
  minHnrowsNxBy2 = wwc.size(0);
  for (k = hnRows; k <= minHnrowsNxBy2; k++) {
    e = wwc[k - 1].re;
    temp_re = fv[k - 1].im;
    twid_re = wwc[k - 1].im;
    twid_im = fv[k - 1].re;
    nd2 = k - hnRows;
    ytmp[nd2].re = e * twid_im + twid_re * temp_re;
    ytmp[nd2].im = e * temp_re - twid_re * twid_im;
  }
  for (i = 0; i < hnRows; i++) {
    double b_ytmp_re_tmp;
    double ytmp_re_tmp;
    minHnrowsNxBy2 = wrapIndex[i];
    e = ytmp[i].re;
    temp_re = reconVar1[i].im;
    twid_re = ytmp[i].im;
    twid_im = reconVar1[i].re;
    temp_im = ytmp[minHnrowsNxBy2 - 1].re;
    b_temp_re_tmp = -ytmp[minHnrowsNxBy2 - 1].im;
    ytmp_re_tmp = reconVar2[i].im;
    b_ytmp_re_tmp = reconVar2[i].re;
    y[i].re = 0.5 * ((e * twid_im - twid_re * temp_re) +
                     (temp_im * b_ytmp_re_tmp - b_temp_re_tmp * ytmp_re_tmp));
    y[i].im = 0.5 * ((e * temp_re + twid_re * twid_im) +
                     (temp_im * ytmp_re_tmp + b_temp_re_tmp * b_ytmp_re_tmp));
    minHnrowsNxBy2 = hnRows + i;
    y[minHnrowsNxBy2].re =
        0.5 * ((e * b_ytmp_re_tmp - twid_re * ytmp_re_tmp) +
               (temp_im * twid_im - b_temp_re_tmp * temp_re));
    y[minHnrowsNxBy2].im =
        0.5 * ((e * ytmp_re_tmp + twid_re * b_ytmp_re_tmp) +
               (temp_im * temp_re + b_temp_re_tmp * twid_im));
  }
}

void FFTImplementationCallback::doHalfLengthRadix2(
    const array<double, 2U> &x, int xoffInit, array<creal_T, 1U> &y,
    int unsigned_nRows, const array<double, 2U> &costab,
    const array<double, 2U> &sintab)
{
  array<creal_T, 1U> reconVar1;
  array<creal_T, 1U> reconVar2;
  array<double, 2U> hcostab;
  array<double, 2U> hsintab;
  array<int, 2U> wrapIndex;
  array<int, 1U> bitrevIndex;
  double b_y_re_tmp;
  double im;
  double re;
  double temp2_im;
  double temp2_re;
  double temp_im;
  double temp_im_tmp;
  double temp_re;
  double temp_re_tmp;
  double y_re_tmp;
  int hszCostab;
  int i;
  int iheight;
  int ihi;
  int istart;
  int j;
  int ju;
  int k;
  int nRowsD2_tmp;
  int nRows_tmp;
  boolean_T tst;
  nRows_tmp = unsigned_nRows / 2;
  iheight = y.size(0);
  if (iheight > nRows_tmp) {
    iheight = nRows_tmp;
  }
  istart = iheight - 2;
  j = nRows_tmp - 2;
  nRowsD2_tmp = nRows_tmp / 2;
  k = nRowsD2_tmp / 2;
  hszCostab = static_cast<int>(static_cast<unsigned int>(costab.size(1)) >> 1);
  hcostab.set_size(1, hszCostab);
  hsintab.set_size(1, hszCostab);
  for (i = 0; i < hszCostab; i++) {
    ju = ((i + 1) << 1) - 2;
    hcostab[i] = costab[ju];
    hsintab[i] = sintab[ju];
  }
  reconVar1.set_size(nRows_tmp);
  reconVar2.set_size(nRows_tmp);
  wrapIndex.set_size(1, nRows_tmp);
  ju = 0;
  hszCostab = 1;
  bitrevIndex.set_size(nRows_tmp);
  for (i = 0; i < nRows_tmp; i++) {
    temp2_re = sintab[i];
    temp2_im = costab[i];
    reconVar1[i].re = temp2_re + 1.0;
    reconVar1[i].im = -temp2_im;
    reconVar2[i].re = 1.0 - temp2_re;
    reconVar2[i].im = temp2_im;
    if (i + 1 != 1) {
      wrapIndex[i] = (nRows_tmp - i) + 1;
    } else {
      wrapIndex[0] = 1;
    }
    bitrevIndex[i] = 0;
  }
  for (ihi = 0; ihi <= istart; ihi++) {
    bitrevIndex[ihi] = hszCostab;
    hszCostab = nRows_tmp;
    tst = true;
    while (tst) {
      hszCostab >>= 1;
      ju ^= hszCostab;
      tst = ((ju & hszCostab) == 0);
    }
    hszCostab = ju + 1;
  }
  bitrevIndex[iheight - 1] = hszCostab;
  if ((static_cast<unsigned int>(x.size(0)) & 1U) == 0U) {
    tst = true;
    iheight = x.size(0);
  } else if (x.size(0) >= unsigned_nRows) {
    tst = true;
    iheight = unsigned_nRows;
  } else {
    tst = false;
    iheight = x.size(0) - 1;
  }
  if (iheight > unsigned_nRows) {
    iheight = unsigned_nRows;
  }
  hszCostab = iheight / 2;
  for (i = 0; i < hszCostab; i++) {
    ju = xoffInit + (i << 1);
    y[bitrevIndex[i] - 1].re = x[ju];
    y[bitrevIndex[i] - 1].im = x[ju + 1];
  }
  if (!tst) {
    if (hszCostab - 1 < 0) {
      iheight = xoffInit;
    } else {
      iheight = xoffInit + (hszCostab << 1);
    }
    y[bitrevIndex[hszCostab] - 1].re = x[iheight];
    y[bitrevIndex[hszCostab] - 1].im = 0.0;
  }
  if (nRows_tmp > 1) {
    for (i = 0; i <= j; i += 2) {
      temp2_re = y[i + 1].re;
      temp2_im = y[i + 1].im;
      temp_re = temp2_re;
      temp_im = temp2_im;
      re = y[i].re;
      im = y[i].im;
      temp2_re = re - temp2_re;
      temp2_im = im - temp2_im;
      y[i + 1].re = temp2_re;
      y[i + 1].im = temp2_im;
      re += temp_re;
      im += temp_im;
      y[i].re = re;
      y[i].im = im;
    }
  }
  hszCostab = 2;
  ju = 4;
  iheight = ((k - 1) << 2) + 1;
  while (k > 0) {
    int b_temp_re_tmp;
    for (i = 0; i < iheight; i += ju) {
      b_temp_re_tmp = i + hszCostab;
      temp_re = y[b_temp_re_tmp].re;
      temp_im = y[b_temp_re_tmp].im;
      y[b_temp_re_tmp].re = y[i].re - temp_re;
      y[b_temp_re_tmp].im = y[i].im - temp_im;
      y[i].re = y[i].re + temp_re;
      y[i].im = y[i].im + temp_im;
    }
    istart = 1;
    for (j = k; j < nRowsD2_tmp; j += k) {
      temp2_re = hcostab[j];
      temp2_im = hsintab[j];
      i = istart;
      ihi = istart + iheight;
      while (i < ihi) {
        b_temp_re_tmp = i + hszCostab;
        temp_re_tmp = y[b_temp_re_tmp].im;
        temp_im = y[b_temp_re_tmp].re;
        temp_re = temp2_re * temp_im - temp2_im * temp_re_tmp;
        temp_im = temp2_re * temp_re_tmp + temp2_im * temp_im;
        y[b_temp_re_tmp].re = y[i].re - temp_re;
        y[b_temp_re_tmp].im = y[i].im - temp_im;
        y[i].re = y[i].re + temp_re;
        y[i].im = y[i].im + temp_im;
        i += ju;
      }
      istart++;
    }
    k /= 2;
    hszCostab = ju;
    ju += ju;
    iheight -= hszCostab;
  }
  temp_re_tmp = y[0].re;
  temp_im = y[0].im;
  y_re_tmp = y[0].re * reconVar1[0].re;
  temp2_im = y[0].re * reconVar1[0].im;
  temp_re = -y[0].im;
  b_y_re_tmp = temp_re_tmp * reconVar2[0].re;
  temp2_re = temp_re_tmp * reconVar2[0].im;
  y[0].re = 0.5 * ((y_re_tmp - y[0].im * reconVar1[0].im) +
                   (b_y_re_tmp - temp_re * reconVar2[0].im));
  y[0].im = 0.5 * ((temp2_im + y[0].im * reconVar1[0].re) +
                   (temp2_re + temp_re * reconVar2[0].re));
  y[nRows_tmp].re = 0.5 * ((b_y_re_tmp - temp_im * reconVar2[0].im) +
                           (y_re_tmp - temp_re * reconVar1[0].im));
  y[nRows_tmp].im = 0.5 * ((temp2_re + temp_im * reconVar2[0].re) +
                           (temp2_im + temp_re * reconVar1[0].re));
  for (i = 2; i <= nRowsD2_tmp; i++) {
    temp_re_tmp = y[i - 1].re;
    temp_im_tmp = y[i - 1].im;
    ju = wrapIndex[i - 1];
    temp2_im = y[ju - 1].re;
    temp_re = y[ju - 1].im;
    y_re_tmp = reconVar1[i - 1].im;
    b_y_re_tmp = reconVar1[i - 1].re;
    temp2_re = reconVar2[i - 1].im;
    temp_im = reconVar2[i - 1].re;
    y[i - 1].re = 0.5 * ((temp_re_tmp * b_y_re_tmp - temp_im_tmp * y_re_tmp) +
                         (temp2_im * temp_im - -temp_re * temp2_re));
    y[i - 1].im = 0.5 * ((temp_re_tmp * y_re_tmp + temp_im_tmp * b_y_re_tmp) +
                         (temp2_im * temp2_re + -temp_re * temp_im));
    hszCostab = (nRows_tmp + i) - 1;
    y[hszCostab].re = 0.5 * ((temp_re_tmp * temp_im - temp_im_tmp * temp2_re) +
                             (temp2_im * b_y_re_tmp - -temp_re * y_re_tmp));
    y[hszCostab].im = 0.5 * ((temp_re_tmp * temp2_re + temp_im_tmp * temp_im) +
                             (temp2_im * y_re_tmp + -temp_re * b_y_re_tmp));
    re = reconVar1[ju - 1].im;
    im = reconVar1[ju - 1].re;
    temp_im = reconVar2[ju - 1].im;
    temp2_re = reconVar2[ju - 1].re;
    y[ju - 1].re = 0.5 * ((temp2_im * im - temp_re * re) +
                          (temp_re_tmp * temp2_re - -temp_im_tmp * temp_im));
    y[ju - 1].im = 0.5 * ((temp2_im * re + temp_re * im) +
                          (temp_re_tmp * temp_im + -temp_im_tmp * temp2_re));
    ju = (ju + nRows_tmp) - 1;
    y[ju].re = 0.5 * ((temp2_im * temp2_re - temp_re * temp_im) +
                      (temp_re_tmp * im - -temp_im_tmp * re));
    y[ju].im = 0.5 * ((temp2_im * temp_im + temp_re * temp2_re) +
                      (temp_re_tmp * re + -temp_im_tmp * im));
  }
  if (nRowsD2_tmp != 0) {
    temp_re_tmp = y[nRowsD2_tmp].re;
    temp_im_tmp = y[nRowsD2_tmp].im;
    y_re_tmp = reconVar1[nRowsD2_tmp].im;
    b_y_re_tmp = reconVar1[nRowsD2_tmp].re;
    temp2_re = temp_re_tmp * b_y_re_tmp;
    temp2_im = temp_re_tmp * y_re_tmp;
    temp_im = reconVar2[nRowsD2_tmp].im;
    re = reconVar2[nRowsD2_tmp].re;
    im = temp_re_tmp * re;
    temp_re = temp_re_tmp * temp_im;
    y[nRowsD2_tmp].re = 0.5 * ((temp2_re - temp_im_tmp * y_re_tmp) +
                               (im - -temp_im_tmp * temp_im));
    y[nRowsD2_tmp].im = 0.5 * ((temp2_im + temp_im_tmp * b_y_re_tmp) +
                               (temp_re + -temp_im_tmp * re));
    ju = nRows_tmp + nRowsD2_tmp;
    y[ju].re = 0.5 * ((im - temp_im_tmp * temp_im) +
                      (temp2_re - -temp_im_tmp * y_re_tmp));
    y[ju].im = 0.5 * ((temp_re + temp_im_tmp * re) +
                      (temp2_im + -temp_im_tmp * b_y_re_tmp));
  }
}

void FFTImplementationCallback::r2br_r2dit_trig_impl(
    const array<creal_T, 1U> &x, int unsigned_nRows,
    const array<double, 2U> &costab, const array<double, 2U> &sintab,
    array<creal_T, 1U> &y)
{
  double im;
  double temp_im;
  double temp_re;
  double temp_re_tmp;
  int i;
  int iDelta2;
  int iheight;
  int iy;
  int ju;
  int k;
  int nRowsD2;
  y.set_size(unsigned_nRows);
  if (unsigned_nRows > x.size(0)) {
    y.set_size(unsigned_nRows);
    for (iy = 0; iy < unsigned_nRows; iy++) {
      y[iy].re = 0.0;
      y[iy].im = 0.0;
    }
  }
  iDelta2 = x.size(0);
  if (iDelta2 > unsigned_nRows) {
    iDelta2 = unsigned_nRows;
  }
  iheight = unsigned_nRows - 2;
  nRowsD2 = unsigned_nRows / 2;
  k = nRowsD2 / 2;
  iy = 0;
  ju = 0;
  for (i = 0; i <= iDelta2 - 2; i++) {
    boolean_T tst;
    y[iy] = x[i];
    iy = unsigned_nRows;
    tst = true;
    while (tst) {
      iy >>= 1;
      ju ^= iy;
      tst = ((ju & iy) == 0);
    }
    iy = ju;
  }
  if (iDelta2 - 2 < 0) {
    iDelta2 = 0;
  } else {
    iDelta2--;
  }
  y[iy] = x[iDelta2];
  if (unsigned_nRows > 1) {
    for (i = 0; i <= iheight; i += 2) {
      temp_re_tmp = y[i + 1].re;
      temp_im = y[i + 1].im;
      temp_re = y[i].re;
      im = y[i].im;
      y[i + 1].re = temp_re - temp_re_tmp;
      y[i + 1].im = y[i].im - y[i + 1].im;
      im += temp_im;
      y[i].re = temp_re + temp_re_tmp;
      y[i].im = im;
    }
  }
  iy = 2;
  iDelta2 = 4;
  iheight = ((k - 1) << 2) + 1;
  while (k > 0) {
    int b_temp_re_tmp;
    for (i = 0; i < iheight; i += iDelta2) {
      b_temp_re_tmp = i + iy;
      temp_re = y[b_temp_re_tmp].re;
      temp_im = y[b_temp_re_tmp].im;
      y[b_temp_re_tmp].re = y[i].re - temp_re;
      y[b_temp_re_tmp].im = y[i].im - temp_im;
      y[i].re = y[i].re + temp_re;
      y[i].im = y[i].im + temp_im;
    }
    ju = 1;
    for (int j{k}; j < nRowsD2; j += k) {
      double twid_im;
      double twid_re;
      int ihi;
      twid_re = costab[j];
      twid_im = sintab[j];
      i = ju;
      ihi = ju + iheight;
      while (i < ihi) {
        b_temp_re_tmp = i + iy;
        temp_re_tmp = y[b_temp_re_tmp].im;
        im = y[b_temp_re_tmp].re;
        temp_re = twid_re * im - twid_im * temp_re_tmp;
        temp_im = twid_re * temp_re_tmp + twid_im * im;
        y[b_temp_re_tmp].re = y[i].re - temp_re;
        y[b_temp_re_tmp].im = y[i].im - temp_im;
        y[i].re = y[i].re + temp_re;
        y[i].im = y[i].im + temp_im;
        i += iDelta2;
      }
      ju++;
    }
    k /= 2;
    iy = iDelta2;
    iDelta2 += iDelta2;
    iheight -= iy;
  }
}

void FFTImplementationCallback::dobluesteinfft(
    const array<double, 2U> &x, int n2blue, int nfft,
    const array<double, 2U> &costab, const array<double, 2U> &sintab,
    const array<double, 2U> &sintabinv, array<creal_T, 2U> &y)
{
  array<creal_T, 1U> fv;
  array<creal_T, 1U> fy;
  array<creal_T, 1U> r;
  array<creal_T, 1U> wwc;
  double temp_im;
  double temp_re;
  double twid_im;
  double twid_re;
  int b_i;
  int b_k;
  int b_y;
  int i;
  int ihi;
  int iy;
  int j;
  int ju;
  int minNrowsNx;
  int nInt2m1;
  int nRowsD2;
  int xoff;
  boolean_T tst;
  if ((nfft != 1) && ((static_cast<unsigned int>(nfft) & 1U) == 0U)) {
    int nInt2;
    int nRows;
    int rt;
    nRows = nfft / 2;
    nInt2m1 = (nRows + nRows) - 1;
    wwc.set_size(nInt2m1);
    rt = 0;
    wwc[nRows - 1].re = 1.0;
    wwc[nRows - 1].im = 0.0;
    nInt2 = nRows << 1;
    for (int k{0}; k <= nRows - 2; k++) {
      double nt_im;
      b_y = ((k + 1) << 1) - 1;
      if (nInt2 - rt <= b_y) {
        rt += b_y - nInt2;
      } else {
        rt += b_y;
      }
      nt_im = -3.1415926535897931 * static_cast<double>(rt) /
              static_cast<double>(nRows);
      i = (nRows - k) - 2;
      wwc[i].re = std::cos(nt_im);
      wwc[i].im = -std::sin(nt_im);
    }
    i = nInt2m1 - 1;
    for (int k{i}; k >= nRows; k--) {
      wwc[k] = wwc[(nInt2m1 - k) - 1];
    }
  } else {
    int nInt2;
    int rt;
    nInt2m1 = (nfft + nfft) - 1;
    wwc.set_size(nInt2m1);
    rt = 0;
    wwc[nfft - 1].re = 1.0;
    wwc[nfft - 1].im = 0.0;
    nInt2 = nfft << 1;
    for (int k{0}; k <= nfft - 2; k++) {
      double nt_im;
      b_y = ((k + 1) << 1) - 1;
      if (nInt2 - rt <= b_y) {
        rt += b_y - nInt2;
      } else {
        rt += b_y;
      }
      nt_im = -3.1415926535897931 * static_cast<double>(rt) /
              static_cast<double>(nfft);
      i = (nfft - k) - 2;
      wwc[i].re = std::cos(nt_im);
      wwc[i].im = -std::sin(nt_im);
    }
    i = nInt2m1 - 1;
    for (int k{i}; k >= nfft; k--) {
      wwc[k] = wwc[(nInt2m1 - k) - 1];
    }
  }
  nInt2m1 = x.size(0);
  y.set_size(nfft, x.size(1));
  if (nfft > x.size(0)) {
    y.set_size(nfft, x.size(1));
    b_y = nfft * x.size(1);
    for (i = 0; i < b_y; i++) {
      y[i].re = 0.0;
      y[i].im = 0.0;
    }
  }
  b_y = x.size(1);
#pragma omp parallel for num_threads(omp_get_max_threads()) private(           \
        fv, fy, r, xoff, minNrowsNx, iy, b_k, j, nRowsD2, ju, b_i, tst,        \
            temp_re, temp_im, twid_re, twid_im, ihi)

  for (int chan = 0; chan < b_y; chan++) {
    xoff = chan * nInt2m1;
    r.set_size(nfft);
    if (nfft > x.size(0)) {
      r.set_size(nfft);
      for (minNrowsNx = 0; minNrowsNx < nfft; minNrowsNx++) {
        r[minNrowsNx].re = 0.0;
        r[minNrowsNx].im = 0.0;
      }
    }
    if ((n2blue != 1) && ((static_cast<unsigned int>(nfft) & 1U) == 0U)) {
      FFTImplementationCallback::doHalfLengthBluestein(
          x, xoff, r, x.size(0), nfft, n2blue, wwc, costab, sintab, costab,
          sintabinv);
    } else {
      minNrowsNx = x.size(0);
      if (nfft <= minNrowsNx) {
        minNrowsNx = nfft;
      }
      for (b_k = 0; b_k < minNrowsNx; b_k++) {
        r[b_k].re = wwc[(nfft + b_k) - 1].re * x[xoff + b_k];
        r[b_k].im = wwc[(nfft + b_k) - 1].im * -x[xoff + b_k];
      }
      minNrowsNx++;
      for (b_k = minNrowsNx; b_k <= nfft; b_k++) {
        r[b_k - 1].re = 0.0;
        r[b_k - 1].im = 0.0;
      }
      fy.set_size(n2blue);
      if (n2blue > r.size(0)) {
        fy.set_size(n2blue);
        for (minNrowsNx = 0; minNrowsNx < n2blue; minNrowsNx++) {
          fy[minNrowsNx].re = 0.0;
          fy[minNrowsNx].im = 0.0;
        }
      }
      iy = r.size(0);
      j = n2blue;
      if (iy <= n2blue) {
        j = iy;
      }
      xoff = n2blue - 2;
      nRowsD2 = n2blue / 2;
      b_k = nRowsD2 / 2;
      iy = 0;
      ju = 0;
      for (b_i = 0; b_i <= j - 2; b_i++) {
        fy[iy] = r[b_i];
        minNrowsNx = n2blue;
        tst = true;
        while (tst) {
          minNrowsNx >>= 1;
          ju ^= minNrowsNx;
          tst = ((ju & minNrowsNx) == 0);
        }
        iy = ju;
      }
      if (j - 2 < 0) {
        minNrowsNx = 0;
      } else {
        minNrowsNx = j - 1;
      }
      fy[iy] = r[minNrowsNx];
      if (n2blue > 1) {
        for (b_i = 0; b_i <= xoff; b_i += 2) {
          temp_re = fy[b_i + 1].re;
          temp_im = fy[b_i + 1].im;
          twid_re = fy[b_i].re;
          twid_im = fy[b_i].im;
          fy[b_i + 1].re = fy[b_i].re - fy[b_i + 1].re;
          fy[b_i + 1].im = fy[b_i].im - fy[b_i + 1].im;
          twid_re += temp_re;
          twid_im += temp_im;
          fy[b_i].re = twid_re;
          fy[b_i].im = twid_im;
        }
      }
      minNrowsNx = 2;
      xoff = 4;
      iy = ((b_k - 1) << 2) + 1;
      while (b_k > 0) {
        for (b_i = 0; b_i < iy; b_i += xoff) {
          temp_re = fy[b_i + minNrowsNx].re;
          temp_im = fy[b_i + minNrowsNx].im;
          fy[b_i + minNrowsNx].re = fy[b_i].re - temp_re;
          fy[b_i + minNrowsNx].im = fy[b_i].im - temp_im;
          fy[b_i].re = fy[b_i].re + temp_re;
          fy[b_i].im = fy[b_i].im + temp_im;
        }
        ju = 1;
        for (j = b_k; j < nRowsD2; j += b_k) {
          twid_re = costab[j];
          twid_im = sintab[j];
          b_i = ju;
          ihi = ju + iy;
          while (b_i < ihi) {
            temp_re = twid_re * fy[b_i + minNrowsNx].re -
                      twid_im * fy[b_i + minNrowsNx].im;
            temp_im = twid_re * fy[b_i + minNrowsNx].im +
                      twid_im * fy[b_i + minNrowsNx].re;
            fy[b_i + minNrowsNx].re = fy[b_i].re - temp_re;
            fy[b_i + minNrowsNx].im = fy[b_i].im - temp_im;
            fy[b_i].re = fy[b_i].re + temp_re;
            fy[b_i].im = fy[b_i].im + temp_im;
            b_i += xoff;
          }
          ju++;
        }
        b_k /= 2;
        minNrowsNx = xoff;
        xoff += xoff;
        iy -= minNrowsNx;
      }
      FFTImplementationCallback::r2br_r2dit_trig_impl(wwc, n2blue, costab,
                                                      sintab, fv);
      iy = fy.size(0);
      for (minNrowsNx = 0; minNrowsNx < iy; minNrowsNx++) {
        twid_im = fy[minNrowsNx].re * fv[minNrowsNx].im +
                  fy[minNrowsNx].im * fv[minNrowsNx].re;
        fy[minNrowsNx].re = fy[minNrowsNx].re * fv[minNrowsNx].re -
                            fy[minNrowsNx].im * fv[minNrowsNx].im;
        fy[minNrowsNx].im = twid_im;
      }
      FFTImplementationCallback::r2br_r2dit_trig_impl(fy, n2blue, costab,
                                                      sintabinv, fv);
      if (fv.size(0) > 1) {
        twid_re = 1.0 / static_cast<double>(fv.size(0));
        iy = fv.size(0);
        for (minNrowsNx = 0; minNrowsNx < iy; minNrowsNx++) {
          fv[minNrowsNx].re = twid_re * fv[minNrowsNx].re;
          fv[minNrowsNx].im = twid_re * fv[minNrowsNx].im;
        }
      }
      minNrowsNx = wwc.size(0);
      for (b_k = nfft; b_k <= minNrowsNx; b_k++) {
        r[b_k - nfft].re =
            wwc[b_k - 1].re * fv[b_k - 1].re + wwc[b_k - 1].im * fv[b_k - 1].im;
        r[b_k - nfft].im =
            wwc[b_k - 1].re * fv[b_k - 1].im - wwc[b_k - 1].im * fv[b_k - 1].re;
      }
    }
    iy = y.size(0);
    for (minNrowsNx = 0; minNrowsNx < iy; minNrowsNx++) {
      y[minNrowsNx + y.size(0) * chan] = r[minNrowsNx];
    }
  }
}

void FFTImplementationCallback::r2br_r2dit_trig(const array<double, 2U> &x,
                                                int n1_unsigned,
                                                const array<double, 2U> &costab,
                                                const array<double, 2U> &sintab,
                                                array<creal_T, 2U> &y)
{
  array<creal_T, 1U> r;
  int i1;
  int loop_ub_tmp;
  int nrows;
  int xoff;
  nrows = x.size(0);
  y.set_size(n1_unsigned, x.size(1));
  if (n1_unsigned > x.size(0)) {
    y.set_size(n1_unsigned, x.size(1));
    loop_ub_tmp = n1_unsigned * x.size(1);
    for (int i{0}; i < loop_ub_tmp; i++) {
      y[i].re = 0.0;
      y[i].im = 0.0;
    }
  }
  loop_ub_tmp = x.size(1);
#pragma omp parallel for num_threads(omp_get_max_threads()) private(r, xoff, i1)

  for (int chan = 0; chan < loop_ub_tmp; chan++) {
    xoff = chan * nrows;
    r.set_size(n1_unsigned);
    if (n1_unsigned > x.size(0)) {
      r.set_size(n1_unsigned);
      for (i1 = 0; i1 < n1_unsigned; i1++) {
        r[i1].re = 0.0;
        r[i1].im = 0.0;
      }
    }
    if (n1_unsigned != 1) {
      FFTImplementationCallback::doHalfLengthRadix2(x, xoff, r, n1_unsigned,
                                                    costab, sintab);
    } else {
      r[0].re = x[xoff];
      r[0].im = 0.0;
    }
    xoff = y.size(0);
    for (i1 = 0; i1 < xoff; i1++) {
      y[i1 + y.size(0) * chan] = r[i1];
    }
  }
}

} // namespace fft
} // namespace internal
} // namespace coder

// End of code generation (FFTImplementationCallback.cpp)
