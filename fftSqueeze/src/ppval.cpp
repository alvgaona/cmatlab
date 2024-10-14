//
// Trial License - for use to evaluate programs for possible purchase as
// an end-user only.
//
// ppval.cpp
//
// Code generation for function 'ppval'
//

// Include files
#include "ppval.h"
#include "fftSqueeze_internal_types.h"
#include "rt_nonfinite.h"
#include "coder_bounded_array.h"
#include "omp.h"
#include <cmath>

// Function Definitions
namespace coder {
int ppval(const struct_T &pp, const double x_data[], double v_data[])
{
  double v;
  double xloc;
  int coefStride;
  int high_i;
  int low_i;
  int low_ip1;
  int mid_i;
  int numTerms;
  int v_size;
  coefStride = pp.breaks.size[1] - 1;
  numTerms = pp.coefs.size[1];
  v_size = 128;
#pragma omp parallel for num_threads(omp_get_max_threads()) private(           \
        high_i, v, low_i, low_ip1, xloc, mid_i)

  for (int ix = 0; ix < 128; ix++) {
    if ((numTerms > 1) && std::isnan(x_data[ix])) {
      v = x_data[ix];
    } else {
      high_i = pp.breaks.size[1];
      low_i = 1;
      low_ip1 = 2;
      while (high_i > low_ip1) {
        mid_i = (low_i >> 1) + (high_i >> 1);
        if (((static_cast<unsigned int>(low_i) & 1U) == 1U) &&
            ((static_cast<unsigned int>(high_i) & 1U) == 1U)) {
          mid_i++;
        }
        if (x_data[ix] >= pp.breaks.data[mid_i - 1]) {
          low_i = mid_i;
          low_ip1 = mid_i + 1;
        } else {
          high_i = mid_i;
        }
      }
      xloc = x_data[ix] - pp.breaks.data[low_i - 1];
      v = pp.coefs.data[low_i - 1];
      for (high_i = 2; high_i <= numTerms; high_i++) {
        v = xloc * v + pp.coefs.data[(low_i + (high_i - 1) * coefStride) - 1];
      }
    }
    v_data[ix] = v;
  }
  return v_size;
}

} // namespace coder

// End of code generation (ppval.cpp)
