//
// Trial License - for use to evaluate programs for possible purchase as
// an end-user only.
//
// imresize.cpp
//
// Code generation for function 'imresize'
//

// Include files
#include "imresize.h"
#include "rt_nonfinite.h"
#include "coder_array.h"
#include "omp.h"
#include <cmath>
#include <emmintrin.h>

// Function Declarations
namespace coder {
static void b_resizeAlongDim2D(const array<unsigned char, 3U> &in,
                               const array<double, 2U> &weights,
                               const array<int, 2U> &indices, double out_length,
                               array<unsigned char, 3U> &out);

static void contributions(int in_length, double out_length,
                          array<double, 2U> &weights, array<int, 2U> &indices);

static void resizeAlongDim2D(const array<unsigned char, 3U> &in,
                             const array<double, 2U> &weights,
                             const array<int, 2U> &indices, double out_length,
                             array<unsigned char, 3U> &out);

} // namespace coder

// Function Definitions
namespace coder {
static void b_resizeAlongDim2D(const array<unsigned char, 3U> &in,
                               const array<double, 2U> &weights,
                               const array<int, 2U> &indices, double out_length,
                               array<unsigned char, 3U> &out)
{
  double sumVal1;
  double v;
  int i;
  int i1;
  int idx;
  int k;
  int linearInds;
  int outCInd;
  int pixelIndex;
  int ub_loop;
  ub_loop = in.size(0);
  for (int pInd{0}; pInd < 3; pInd++) {
    int b_pInd;
    b_pInd = pInd + 1;
#pragma omp parallel for num_threads(omp_get_max_threads()) private(           \
        pixelIndex, linearInds, sumVal1, idx, i, outCInd, i1, k, v)

    for (int inRInd = 0; inRInd < ub_loop; inRInd++) {
      idx = (inRInd + in.size(0) * in.size(1) * (b_pInd - 1)) + 1;
      i = static_cast<int>(out_length);
      for (outCInd = 0; outCInd < i; outCInd++) {
        sumVal1 = 0.0;
        //  Core - second dimension
        i1 = weights.size(0);
        linearInds = weights.size(0) * outCInd + 1;
        for (k = 0; k < i1; k++) {
          pixelIndex = idx + (indices[(linearInds + k) - 1] - 1) * in.size(0);
          sumVal1 += weights[(linearInds + k) - 1] *
                     static_cast<double>(in[pixelIndex - 1]);
        }
        v = std::abs(sumVal1);
        if (v < 4.503599627370496E+15) {
          if (v >= 0.5) {
            v = std::floor(sumVal1 + 0.5);
          } else {
            v = sumVal1 * 0.0;
          }
        } else {
          v = sumVal1;
        }
        if (sumVal1 > 255.0) {
          out[(inRInd + out.size(0) * outCInd) +
              out.size(0) * out.size(1) * (b_pInd - 1)] = MAX_uint8_T;
        } else if (sumVal1 < 0.0) {
          out[(inRInd + out.size(0) * outCInd) +
              out.size(0) * out.size(1) * (b_pInd - 1)] = 0U;
        } else {
          out[(inRInd + out.size(0) * outCInd) +
              out.size(0) * out.size(1) * (b_pInd - 1)] =
              static_cast<unsigned char>(v);
        }
      }
    }
  }
}

static void contributions(int in_length, double out_length,
                          array<double, 2U> &weights, array<int, 2U> &indices)
{
  __m128d r;
  __m128d r1;
  array<double, 2U> absx;
  array<double, 2U> absx2;
  array<double, 2U> absx3;
  array<double, 2U> b_x;
  array<double, 2U> y;
  array<double, 1U> u;
  array<double, 1U> x;
  array<int, 2U> aux;
  array<int, 2U> b_indices;
  array<int, 1U> left;
  int acoef;
  int b_k;
  int bcoef;
  int csz_idx_0_tmp;
  int i;
  int loop_ub_tmp;
  int nx_tmp;
  int xoffset;
  signed char tmp_data[10];
  boolean_T copyCols[10];
  //  Contributions, using pixel indices
  if (out_length < 1.0) {
    y.set_size(1, 0);
  } else {
    y.set_size(1, static_cast<int>(out_length - 1.0) + 1);
    xoffset = static_cast<int>(out_length - 1.0);
    for (i = 0; i <= xoffset; i++) {
      y[i] = static_cast<double>(i) + 1.0;
    }
  }
  xoffset = y.size(1);
  u.set_size(y.size(1));
  acoef = (y.size(1) / 2) << 1;
  bcoef = acoef - 2;
  for (i = 0; i <= bcoef; i += 2) {
    r = _mm_loadu_pd(&y[i]);
    r1 = _mm_set1_pd(0.5);
    _mm_storeu_pd(&u[i], _mm_sub_pd(_mm_div_pd(r, r1), r1));
  }
  for (i = acoef; i < xoffset; i++) {
    u[i] = y[i] / 0.5 - 0.5;
  }
  xoffset = u.size(0);
  x.set_size(u.size(0));
  acoef = (u.size(0) / 2) << 1;
  bcoef = acoef - 2;
  for (i = 0; i <= bcoef; i += 2) {
    r = _mm_loadu_pd(&u[i]);
    _mm_storeu_pd(&x[i], _mm_sub_pd(r, _mm_set1_pd(4.0)));
  }
  for (i = acoef; i < xoffset; i++) {
    x[i] = u[i] - 4.0;
  }
  nx_tmp = x.size(0);
  for (int k{0}; k < nx_tmp; k++) {
    x[k] = std::floor(x[k]);
  }
  left.set_size(x.size(0));
  for (i = 0; i < nx_tmp; i++) {
    left[i] = static_cast<int>(x[i]);
  }
  b_indices.set_size(x.size(0), 10);
  if (left.size(0) != 0) {
    acoef = (left.size(0) != 1);
    for (int k{0}; k < 10; k++) {
      for (b_k = 0; b_k < nx_tmp; b_k++) {
        b_indices[b_k + b_indices.size(0) * k] = left[acoef * b_k] + k;
      }
    }
  }
  absx.set_size(x.size(0), 10);
  loop_ub_tmp = b_indices.size(0) * 10;
  for (i = 0; i < loop_ub_tmp; i++) {
    absx[i] = b_indices[i];
  }
  xoffset = absx.size(0);
  acoef = u.size(0);
  if (xoffset <= acoef) {
    acoef = xoffset;
  }
  if (absx.size(0) == 1) {
    xoffset = u.size(0);
  } else if (u.size(0) == 1) {
    xoffset = absx.size(0);
  } else if (u.size(0) == absx.size(0)) {
    xoffset = u.size(0);
  } else {
    xoffset = acoef;
  }
  b_x.set_size(xoffset, 10);
  if (xoffset != 0) {
    acoef = (u.size(0) != 1);
    bcoef = (absx.size(0) != 1);
    for (int k{0}; k < 10; k++) {
      for (b_k = 0; b_k < xoffset; b_k++) {
        b_x[b_k + b_x.size(0) * k] =
            u[acoef * b_k] - absx[bcoef * b_k + absx.size(0) * k];
      }
    }
  }
  xoffset = b_x.size(0) * 10;
  b_x.set_size(b_x.size(0), 10);
  acoef = (xoffset / 2) << 1;
  bcoef = acoef - 2;
  for (i = 0; i <= bcoef; i += 2) {
    r = _mm_loadu_pd(&b_x[i]);
    _mm_storeu_pd(&b_x[i], _mm_mul_pd(_mm_set1_pd(0.5), r));
  }
  for (i = acoef; i < xoffset; i++) {
    b_x[i] = 0.5 * b_x[i];
  }
  b_k = b_x.size(0) * 10;
  i = b_x.size(0);
  absx.set_size(b_x.size(0), 10);
  for (int k{0}; k < b_k; k++) {
    absx[k] = std::abs(b_x[k]);
  }
  absx2.set_size(b_x.size(0), 10);
  acoef = (b_k / 2) << 1;
  bcoef = acoef - 2;
  for (int k{0}; k <= bcoef; k += 2) {
    r = _mm_loadu_pd(&absx[k]);
    r1 = _mm_loadu_pd(&absx[k]);
    _mm_storeu_pd(&absx2[k], _mm_mul_pd(r, r1));
  }
  for (int k{acoef}; k < b_k; k++) {
    absx2[k] = absx[k] * absx[k];
  }
  absx3.set_size(b_x.size(0), 10);
  for (int k{0}; k < b_k; k++) {
    if (std::isnan(absx[k])) {
      absx3[k] = rtNaN;
    } else {
      absx3[k] = std::pow(absx[k], 3.0);
    }
  }
  xoffset = absx2.size(0) * 10;
  absx2.set_size(absx2.size(0), 10);
  acoef = (xoffset / 2) << 1;
  bcoef = acoef - 2;
  for (int xj{0}; xj <= bcoef; xj += 2) {
    r = _mm_loadu_pd(&absx2[xj]);
    _mm_storeu_pd(&absx2[xj], _mm_mul_pd(_mm_set1_pd(2.5), r));
  }
  for (int xj{acoef}; xj < xoffset; xj++) {
    absx2[xj] = 2.5 * absx2[xj];
  }
  absx3.set_size(absx3.size(0), 10);
  for (int xj{0}; xj < b_k; xj++) {
    absx3[xj] =
        0.5 * (((1.5 * absx3[xj] - absx2[xj]) + 1.0) *
                   static_cast<double>(absx[xj] <= 1.0) +
               (((-0.5 * absx3[xj] + absx2[xj]) - 4.0 * absx[xj]) + 2.0) *
                   static_cast<double>((absx[xj] > 1.0) && (absx[xj] <= 2.0)));
  }
  if (absx3.size(0) == 0) {
    u.set_size(0);
  } else {
    u.set_size(b_x.size(0));
    for (int xj{0}; xj < i; xj++) {
      u[xj] = absx3[xj];
    }
    acoef = (i / 2) << 1;
    bcoef = acoef - 2;
    for (int k{0}; k < 9; k++) {
      xoffset = (k + 1) * i;
      for (int xj{0}; xj <= bcoef; xj += 2) {
        r = _mm_loadu_pd(&u[xj]);
        r1 = _mm_loadu_pd(&absx3[xoffset + xj]);
        _mm_storeu_pd(&u[xj], _mm_add_pd(r, r1));
      }
      for (int xj{acoef}; xj < i; xj++) {
        u[xj] = u[xj] + absx3[xoffset + xj];
      }
    }
  }
  absx.set_size(b_x.size(0), 10);
  for (i = 0; i < b_k; i++) {
    absx[i] = absx3[i];
  }
  xoffset = u.size(0);
  acoef = absx3.size(0);
  if (xoffset <= acoef) {
    acoef = xoffset;
  }
  if (u.size(0) == 1) {
    csz_idx_0_tmp = absx3.size(0);
  } else if (absx3.size(0) == 1) {
    csz_idx_0_tmp = u.size(0);
  } else if (absx3.size(0) == u.size(0)) {
    csz_idx_0_tmp = absx3.size(0);
  } else {
    csz_idx_0_tmp = acoef;
  }
  absx3.set_size(csz_idx_0_tmp, 10);
  if (csz_idx_0_tmp != 0) {
    acoef = (absx.size(0) != 1);
    bcoef = (u.size(0) != 1);
    for (int k{0}; k < 10; k++) {
      for (b_k = 0; b_k < csz_idx_0_tmp; b_k++) {
        absx3[b_k + absx3.size(0) * k] =
            absx[acoef * b_k + absx.size(0) * k] / u[bcoef * b_k];
      }
    }
  }
  //  Create the auxiliary matrix:
  xoffset = in_length << 1;
  aux.set_size(1, xoffset);
  aux[0] = 1;
  aux[in_length] = in_length;
  for (b_k = 2; b_k <= in_length; b_k++) {
    aux[b_k - 1] = aux[b_k - 2] + 1;
    acoef = in_length + b_k;
    aux[acoef - 1] = aux[acoef - 2] - 1;
  }
  //  Mirror the out-of-bounds indices using mod:
  for (b_k = 0; b_k < loop_ub_tmp; b_k++) {
    double c_k;
    c_k = static_cast<double>(b_indices[b_k]) - 1.0;
    if (xoffset == 0) {
      if (c_k == 0.0) {
        c_k = 0.0;
      }
    } else if (c_k == 0.0) {
      c_k = 0.0;
    } else {
      c_k = std::fmod(c_k, static_cast<double>(xoffset));
      if (c_k == 0.0) {
        c_k = 0.0;
      } else if (c_k < 0.0) {
        c_k += static_cast<double>(xoffset);
      }
    }
    b_indices[b_k] = aux[static_cast<int>(c_k)];
  }
  for (i = 0; i < 10; i++) {
    copyCols[i] = false;
  }
  acoef = 0;
  for (b_k = 0; b_k < 10; b_k++) {
    boolean_T exitg1;
    xoffset = acoef + absx3.size(0);
    bcoef = acoef;
    acoef += absx3.size(0);
    exitg1 = false;
    while ((!exitg1) && (bcoef + 1 <= xoffset)) {
      if ((absx3[bcoef] == 0.0) || std::isnan(absx3[bcoef])) {
        bcoef++;
      } else {
        copyCols[b_k] = true;
        exitg1 = true;
      }
    }
  }
  acoef = 0;
  xoffset = 0;
  for (b_k = 0; b_k < 10; b_k++) {
    if (copyCols[b_k]) {
      acoef++;
      tmp_data[xoffset] = static_cast<signed char>(b_k);
      xoffset++;
    }
  }
  weights.set_size(acoef, csz_idx_0_tmp);
  for (i = 0; i < csz_idx_0_tmp; i++) {
    for (int xj{0}; xj < acoef; xj++) {
      weights[xj + weights.size(0) * i] =
          absx3[i + absx3.size(0) * tmp_data[xj]];
    }
  }
  indices.set_size(acoef, x.size(0));
  for (i = 0; i < nx_tmp; i++) {
    for (int xj{0}; xj < acoef; xj++) {
      indices[xj + indices.size(0) * i] =
          b_indices[i + b_indices.size(0) * tmp_data[xj]];
    }
  }
}

static void resizeAlongDim2D(const array<unsigned char, 3U> &in,
                             const array<double, 2U> &weights,
                             const array<int, 2U> &indices, double out_length,
                             array<unsigned char, 3U> &out)
{
  double sumVal1;
  int i;
  int i1;
  int in_idx_0;
  int k;
  int linearInds;
  int outRInd;
  int ub_loop;
  ub_loop = static_cast<int>(static_cast<double>(in.size(0) * in.size(1) * 3) /
                             static_cast<double>(in.size(0)));
#pragma omp parallel for num_threads(omp_get_max_threads()) private(           \
        linearInds, sumVal1, i, outRInd, i1, in_idx_0, k)

  for (int inCInd = 0; inCInd < ub_loop; inCInd++) {
    i = static_cast<int>(out_length);
    for (outRInd = 0; outRInd < i; outRInd++) {
      sumVal1 = 0.0;
      i1 = weights.size(0);
      linearInds = weights.size(0) * outRInd + 1;
      //  Core - first dimension
      if (weights.size(0) - 1 >= 0) {
        in_idx_0 = in.size(0);
      }
      for (k = 0; k < i1; k++) {
        sumVal1 +=
            weights[(linearInds + k) - 1] *
            static_cast<double>(
                in[(indices[(linearInds + k) - 1] + in_idx_0 * inCInd) - 1]);
      }
      if (sumVal1 > 255.0) {
        in_idx_0 = out.size(0);
        out[outRInd + in_idx_0 * inCInd] = MAX_uint8_T;
      } else if (sumVal1 < 0.0) {
        in_idx_0 = out.size(0);
        out[outRInd + in_idx_0 * inCInd] = 0U;
      } else {
        in_idx_0 = out.size(0);
        if (sumVal1 < 4.503599627370496E+15) {
          if (sumVal1 >= 0.5) {
            sumVal1 = std::floor(sumVal1 + 0.5);
          } else {
            sumVal1 *= 0.0;
          }
        }
        out[outRInd + in_idx_0 * inCInd] = static_cast<unsigned char>(sumVal1);
      }
    }
  }
}

void imresize(const array<unsigned char, 3U> &Ain,
              array<unsigned char, 3U> &Bout)
{
  array<double, 2U> weights;
  array<int, 2U> indices;
  array<unsigned char, 3U> out;
  //  Resize first dimension
  contributions(Ain.size(0), std::ceil(static_cast<double>(Ain.size(0)) * 0.5),
                weights, indices);
  out.set_size(weights.size(1), Ain.size(1), 3);
  resizeAlongDim2D(Ain, weights, indices, static_cast<double>(weights.size(1)),
                   out);
  //  Resize second dimension
  contributions(Ain.size(1), std::ceil(static_cast<double>(Ain.size(1)) * 0.5),
                weights, indices);
  Bout.set_size(out.size(0), weights.size(1), 3);
  b_resizeAlongDim2D(out, weights, indices,
                     static_cast<double>(weights.size(1)), Bout);
}

} // namespace coder

// End of code generation (imresize.cpp)
