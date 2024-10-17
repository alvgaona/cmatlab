//
// Trial License - for use to evaluate programs for possible purchase as
// an end-user only.
//
// dtwin.cpp
//
// Code generation for function 'dtwin'
//

// Include files
#include "dtwin.h"
#include "fftsqueeze_internal_types.h"
#include "ppval.h"
#include "rt_nonfinite.h"
#include "coder_array.h"
#include <cstring>
#include <emmintrin.h>

// Function Definitions
namespace coder {
namespace b_signal {
namespace internal {
namespace spectral {
void dtwin(const array<double, 1U> &w, double Fs, array<double, 1U> &Wdt)
{
  array<double, 2U> dvdf;
  array<double, 2U> md;
  array<double, 2U> pp_coefs;
  array<double, 2U> s;
  array<double, 1U> b_y;
  array<int, 2U> dx;
  array<unsigned int, 2U> t0_breaks;
  array<unsigned int, 2U> y;
  struct_T expl_temp;
  double xv_data[4];
  double szdvdf_idx_1;
  int b_loop_ub;
  int d31;
  int dnnm2;
  int loop_ub;
  int stride;
  int szs_idx_1;
  signed char outsize_idx_0;
  boolean_T has_endslopes;
  if (w.size(0) < 1) {
    y.set_size(1, 0);
  } else {
    y.set_size(1, w.size(0));
    loop_ub = w.size(0) - 1;
    for (d31 = 0; d31 <= loop_ub; d31++) {
      y[d31] = static_cast<unsigned int>(d31) + 1U;
    }
  }
  has_endslopes = (w.size(0) == y.size(1) + 2);
  if ((y.size(1) <= 2) || ((y.size(1) <= 3) && (!has_endslopes))) {
    has_endslopes = (w.size(0) == y.size(1) + 2);
    if (y.size(1) <= 2) {
      if (has_endslopes) {
        stride = 4;
      } else {
        stride = 2;
      }
    } else {
      stride = 3;
    }
    pp_coefs.set_size(1, stride);
    if (y.size(1) <= 2) {
      if (has_endslopes) {
        if (y.size(1) - 2 >= 0) {
          double dzzdx;
          szs_idx_1 = static_cast<int>(y[1]) - static_cast<int>(y[0]);
          szdvdf_idx_1 = (w[2] - w[1]) / static_cast<double>(szs_idx_1);
          dzzdx = (szdvdf_idx_1 - w[0]) / static_cast<double>(szs_idx_1);
          szdvdf_idx_1 = (w[w.size(0) - 1] - szdvdf_idx_1) /
                         static_cast<double>(szs_idx_1);
          xv_data[0] = (szdvdf_idx_1 - dzzdx) / static_cast<double>(szs_idx_1);
          xv_data[1] = 2.0 * dzzdx - szdvdf_idx_1;
          xv_data[2] = w[0];
          xv_data[3] = w[1];
        }
        pp_coefs.set_size(1, stride);
        for (d31 = 0; d31 < stride; d31++) {
          pp_coefs[d31] = xv_data[d31];
        }
      } else {
        pp_coefs[0] =
            (w[1] - w[0]) / static_cast<double>(static_cast<int>(y[1]) -
                                                static_cast<int>(y[0]));
        pp_coefs[1] = w[0];
      }
      loop_ub = y.size(1);
      t0_breaks.set_size(1, y.size(1));
      for (d31 = 0; d31 < loop_ub; d31++) {
        t0_breaks[d31] = y[d31];
      }
    } else {
      szs_idx_1 = static_cast<int>(y[1]) - static_cast<int>(y[0]);
      szdvdf_idx_1 = (w[1] - w[0]) / static_cast<double>(szs_idx_1);
      pp_coefs[0] =
          ((w[2] - w[1]) / static_cast<double>(static_cast<int>(y[2]) -
                                               static_cast<int>(y[1])) -
           szdvdf_idx_1) /
          static_cast<double>(static_cast<int>(y[2]) - static_cast<int>(y[0]));
      pp_coefs[1] = szdvdf_idx_1 - pp_coefs[0] * static_cast<double>(szs_idx_1);
      pp_coefs[2] = w[0];
      t0_breaks.set_size(1, 2);
      t0_breaks[0] = y[0];
      t0_breaks[1] = y[2];
    }
  } else {
    int nxm1_tmp;
    int yoffset;
    nxm1_tmp = y.size(1) - 1;
    if (has_endslopes) {
      szdvdf_idx_1 = static_cast<double>(w.size(0)) - 3.0;
      szs_idx_1 = w.size(0) - 2;
      yoffset = 1;
    } else {
      szdvdf_idx_1 = static_cast<double>(w.size(0)) - 1.0;
      szs_idx_1 = w.size(0);
      yoffset = 0;
    }
    s.set_size(1, szs_idx_1);
    dvdf.set_size(1, static_cast<int>(szdvdf_idx_1));
    dx.set_size(1, y.size(1) - 1);
    for (int k{0}; k < nxm1_tmp; k++) {
      d31 = static_cast<int>(y[k + 1]) - static_cast<int>(y[k]);
      dx[k] = d31;
      szs_idx_1 = yoffset + k;
      dvdf[k] = (w[szs_idx_1 + 1] - w[szs_idx_1]) / static_cast<double>(d31);
    }
    for (int k{2}; k <= nxm1_tmp; k++) {
      s[k - 1] = 3.0 * (static_cast<double>(dx[k - 1]) * dvdf[k - 2] +
                        static_cast<double>(dx[k - 2]) * dvdf[k - 1]);
    }
    if (has_endslopes) {
      d31 = 0;
      dnnm2 = 0;
      s[0] = w[0] * static_cast<double>(dx[1]);
      s[y.size(1) - 1] =
          static_cast<double>(dx[y.size(1) - 3]) * w[y.size(1) + 1];
    } else {
      d31 = static_cast<int>(y[2]) - static_cast<int>(y[0]);
      dnnm2 = static_cast<int>(y[y.size(1) - 1]) -
              static_cast<int>(y[y.size(1) - 3]);
      szs_idx_1 = dx[0];
      s[0] =
          ((static_cast<double>(szs_idx_1) + 2.0 * static_cast<double>(d31)) *
               static_cast<double>(dx[1]) * dvdf[0] +
           static_cast<double>(szs_idx_1) * static_cast<double>(szs_idx_1) *
               dvdf[1]) /
          static_cast<double>(d31);
      szs_idx_1 = dx[y.size(1) - 2];
      s[y.size(1) - 1] =
          ((static_cast<double>(szs_idx_1) + 2.0 * static_cast<double>(dnnm2)) *
               static_cast<double>(dx[y.size(1) - 3]) * dvdf[y.size(1) - 2] +
           static_cast<double>(szs_idx_1) * static_cast<double>(szs_idx_1) *
               dvdf[y.size(1) - 3]) /
          static_cast<double>(dnnm2);
    }
    loop_ub = y.size(1);
    md.set_size(1, y.size(1));
    szs_idx_1 = dx[1];
    md[0] = szs_idx_1;
    stride = dx[y.size(1) - 3];
    md[y.size(1) - 1] = stride;
    for (int k{2}; k <= nxm1_tmp; k++) {
      md[k - 1] = 2.0 * (static_cast<double>(dx[k - 1]) +
                         static_cast<double>(dx[k - 2]));
    }
    szdvdf_idx_1 = static_cast<double>(szs_idx_1) / md[0];
    md[1] = md[1] - szdvdf_idx_1 * static_cast<double>(d31);
    s[1] = s[1] - szdvdf_idx_1 * s[0];
    for (int k{3}; k <= nxm1_tmp; k++) {
      szdvdf_idx_1 = static_cast<double>(dx[k - 1]) / md[k - 2];
      md[k - 1] = md[k - 1] - szdvdf_idx_1 * static_cast<double>(dx[k - 3]);
      s[k - 1] = s[k - 1] - szdvdf_idx_1 * s[k - 2];
    }
    szdvdf_idx_1 = static_cast<double>(dnnm2) / md[y.size(1) - 2];
    md[y.size(1) - 1] =
        md[y.size(1) - 1] - szdvdf_idx_1 * static_cast<double>(stride);
    s[y.size(1) - 1] = s[y.size(1) - 1] - szdvdf_idx_1 * s[y.size(1) - 2];
    s[y.size(1) - 1] = s[y.size(1) - 1] / md[y.size(1) - 1];
    for (int k{nxm1_tmp}; k >= 2; k--) {
      s[k - 1] = (s[k - 1] - static_cast<double>(dx[k - 2]) * s[k]) / md[k - 1];
    }
    s[0] = (s[0] - static_cast<double>(d31) * s[1]) / md[0];
    szs_idx_1 = s.size(1) - 1;
    pp_coefs.set_size(s.size(1) - 1, 4);
    for (dnnm2 = 0; dnnm2 <= loop_ub - 2; dnnm2++) {
      double d;
      double dzzdx;
      szdvdf_idx_1 = dvdf[dnnm2];
      d = s[dnnm2];
      d31 = dx[dnnm2];
      dzzdx = (szdvdf_idx_1 - d) / static_cast<double>(d31);
      szdvdf_idx_1 = (s[dnnm2 + 1] - szdvdf_idx_1) / static_cast<double>(d31);
      pp_coefs[dnnm2] = (szdvdf_idx_1 - dzzdx) / static_cast<double>(d31);
      pp_coefs[szs_idx_1 + dnnm2] = 2.0 * dzzdx - szdvdf_idx_1;
      pp_coefs[(szs_idx_1 << 1) + dnnm2] = d;
      pp_coefs[3 * szs_idx_1 + dnnm2] = w[yoffset + dnnm2];
    }
    t0_breaks.set_size(1, y.size(1));
    for (d31 = 0; d31 < loop_ub; d31++) {
      t0_breaks[d31] = y[d31];
    }
  }
  expl_temp.coefs.set_size(pp_coefs.size(0), pp_coefs.size(1) - 1);
  szs_idx_1 = pp_coefs.size(0) * (pp_coefs.size(1) - 1);
  for (d31 = 0; d31 < szs_idx_1; d31++) {
    expl_temp.coefs[d31] = 0.0;
  }
  szs_idx_1 = pp_coefs.size(1);
  stride = pp_coefs.size(0);
  if (pp_coefs.size(0) - 1 >= 0) {
    outsize_idx_0 = static_cast<signed char>(pp_coefs.size(1));
    b_loop_ub = pp_coefs.size(1);
  }
  for (dnnm2 = 0; dnnm2 < stride; dnnm2++) {
    d31 = outsize_idx_0;
    std::memset(&xv_data[0], 0,
                static_cast<unsigned int>(b_loop_ub) * sizeof(double));
    for (int k{0}; k < szs_idx_1; k++) {
      xv_data[k] = pp_coefs[dnnm2 + k * stride];
    }
    for (int k{0}; k <= d31 - 2; k++) {
      expl_temp.coefs[dnnm2 + k * stride] =
          xv_data[k] * (static_cast<double>(outsize_idx_0 - k) - 1.0);
    }
  }
  if (w.size(0) < 1) {
    y.set_size(1, 0);
  } else {
    y.set_size(1, w.size(0));
    loop_ub = w.size(0) - 1;
    for (d31 = 0; d31 <= loop_ub; d31++) {
      y[d31] = static_cast<unsigned int>(d31) + 1U;
    }
  }
  loop_ub = t0_breaks.size(1);
  expl_temp.breaks.set_size(1, t0_breaks.size(1));
  for (d31 = 0; d31 < loop_ub; d31++) {
    expl_temp.breaks[d31] = t0_breaks[d31];
  }
  loop_ub = y.size(1);
  b_y.set_size(y.size(1));
  for (d31 = 0; d31 < loop_ub; d31++) {
    b_y[d31] = y[d31];
  }
  ppval(expl_temp, b_y, Wdt);
  szdvdf_idx_1 = Fs / 6.2831853071795862;
  loop_ub = Wdt.size(0);
  szs_idx_1 = (Wdt.size(0) / 2) << 1;
  stride = szs_idx_1 - 2;
  for (d31 = 0; d31 <= stride; d31 += 2) {
    __m128d r;
    r = _mm_loadu_pd(&Wdt[d31]);
    _mm_storeu_pd(&Wdt[d31], _mm_mul_pd(r, _mm_set1_pd(szdvdf_idx_1)));
  }
  for (d31 = szs_idx_1; d31 < loop_ub; d31++) {
    Wdt[d31] = Wdt[d31] * szdvdf_idx_1;
  }
}

} // namespace spectral
} // namespace internal
} // namespace b_signal
} // namespace coder

// End of code generation (dtwin.cpp)
