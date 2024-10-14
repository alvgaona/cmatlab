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
#include "fftSqueeze_internal_types.h"
#include "ppval.h"
#include "rt_nonfinite.h"
#include "coder_bounded_array.h"
#include <cstring>
#include <emmintrin.h>

// Function Definitions
namespace coder {
namespace b_signal {
namespace internal {
namespace spectral {
int dtwin(const double w_data[], double Fs, double Wdt_data[])
{
  __m128d r;
  __m128d r1;
  struct_T expl_temp;
  double pp_coefs_data[508];
  double md_data[128];
  double s_data[128];
  double dvdf_data[127];
  double b_r;
  int Wdt_size;
  int j;
  int k;
  unsigned char pp_breaks_data[128];
  for (k = 0; k <= 124; k += 2) {
    r = _mm_loadu_pd(&w_data[k + 1]);
    r1 = _mm_loadu_pd(&w_data[k]);
    r = _mm_sub_pd(r, r1);
    _mm_storeu_pd(&dvdf_data[k], r);
  }
  dvdf_data[126] = w_data[127] - w_data[126];
  s_data[0] = (5.0 * dvdf_data[0] + dvdf_data[1]) / 2.0;
  s_data[127] = (5.0 * dvdf_data[126] + dvdf_data[125]) / 2.0;
  md_data[0] = 1.0;
  md_data[127] = 1.0;
  for (k = 0; k <= 124; k += 2) {
    r = _mm_loadu_pd(&dvdf_data[k]);
    r1 = _mm_loadu_pd(&dvdf_data[k + 1]);
    r = _mm_add_pd(r, r1);
    r = _mm_mul_pd(_mm_set1_pd(3.0), r);
    _mm_storeu_pd(&s_data[k + 1], r);
    _mm_storeu_pd(&md_data[k + 1], _mm_set1_pd(4.0));
  }
  b_r = 1.0 / md_data[0];
  md_data[1] -= b_r * 2.0;
  s_data[1] -= b_r * s_data[0];
  for (k = 0; k < 125; k++) {
    b_r = 1.0 / md_data[k + 1];
    md_data[k + 2] -= b_r;
    s_data[k + 2] -= b_r * s_data[k + 1];
  }
  b_r = 2.0 / md_data[126];
  md_data[127] -= b_r;
  s_data[127] -= b_r * s_data[126];
  s_data[127] /= md_data[127];
  for (k = 125; k >= 0; k--) {
    s_data[k + 1] = (s_data[k + 1] - s_data[k + 2]) / md_data[k + 1];
  }
  s_data[0] = (s_data[0] - 2.0 * s_data[1]) / md_data[0];
  for (j = 0; j <= 124; j += 2) {
    __m128d r2;
    __m128d r3;
    r = _mm_loadu_pd(&dvdf_data[j]);
    r1 = _mm_loadu_pd(&s_data[j]);
    r2 = _mm_sub_pd(r, r1);
    r3 = _mm_loadu_pd(&s_data[j + 1]);
    r = _mm_sub_pd(r3, r);
    r3 = _mm_sub_pd(r, r2);
    _mm_storeu_pd(&pp_coefs_data[j], r3);
    r2 = _mm_mul_pd(_mm_set1_pd(2.0), r2);
    r = _mm_sub_pd(r2, r);
    _mm_storeu_pd(&pp_coefs_data[j + 127], r);
    _mm_storeu_pd(&pp_coefs_data[j + 254], r1);
    r = _mm_loadu_pd(&w_data[j]);
    _mm_storeu_pd(&pp_coefs_data[j + 381], r);
  }
  double d;
  double dzzdx;
  b_r = dvdf_data[126];
  d = s_data[126];
  dzzdx = b_r - d;
  b_r = s_data[127] - b_r;
  pp_coefs_data[126] = b_r - dzzdx;
  pp_coefs_data[253] = 2.0 * dzzdx - b_r;
  pp_coefs_data[380] = d;
  pp_coefs_data[507] = w_data[126];
  for (int i{0}; i < 128; i++) {
    pp_breaks_data[i] =
        static_cast<unsigned char>(static_cast<unsigned int>(i) + 1U);
  }
  expl_temp.coefs.size[0] = 127;
  expl_temp.coefs.size[1] = 3;
  std::memset(&expl_temp.coefs.data[0], 0, 381U * sizeof(double));
  for (j = 0; j < 127; j++) {
    double xv_data[4];
    for (k = 0; k < 4; k++) {
      xv_data[k] = pp_coefs_data[j + k * 127];
    }
    for (k = 0; k < 3; k++) {
      expl_temp.coefs.data[j + k * 127] =
          xv_data[k] * (3.0 - static_cast<double>(k));
    }
  }
  expl_temp.breaks.size[0] = 1;
  expl_temp.breaks.size[1] = 128;
  for (int i{0}; i < 128; i++) {
    expl_temp.breaks.data[i] = pp_breaks_data[i];
    s_data[i] = static_cast<double>(i) + 1.0;
  }
  Wdt_size = ppval(expl_temp, s_data, Wdt_data);
  b_r = Fs / 6.2831853071795862;
  k = (Wdt_size / 2) << 1;
  j = k - 2;
  for (int i{0}; i <= j; i += 2) {
    r = _mm_loadu_pd(&Wdt_data[i]);
    _mm_storeu_pd(&Wdt_data[i], _mm_mul_pd(r, _mm_set1_pd(b_r)));
  }
  for (int i{k}; i < Wdt_size; i++) {
    Wdt_data[i] *= b_r;
  }
  return Wdt_size;
}

} // namespace spectral
} // namespace internal
} // namespace b_signal
} // namespace coder

// End of code generation (dtwin.cpp)
