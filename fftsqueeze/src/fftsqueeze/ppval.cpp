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
#include "fftsqueeze_internal_types.h"
#include "rt_nonfinite.h"
#include "coder_array.h"
#include "omp.h"
#include <cmath>

// Function Definitions
namespace coder {
void ppval(const struct_T &pp, const array<double, 1U> &x, array<double, 1U> &v)
{
  double b_v;
  double xloc;
  int coefStride;
  int high_i;
  int low_i;
  int low_ip1;
  int mid_i;
  int numTerms;
  int ub_loop;
  coefStride = pp.breaks.size(1) - 1;
  numTerms = pp.coefs.size(1);
  v.set_size(x.size(0));
  ub_loop = x.size(0);
#pragma omp parallel for num_threads(omp_get_max_threads()) private(           \
        high_i, b_v, low_i, low_ip1, xloc, mid_i)

  for (int ix = 0; ix < ub_loop; ix++) {
    if ((numTerms > 1) && std::isnan(x[ix])) {
      b_v = x[ix];
    } else {
      high_i = pp.breaks.size(1);
      low_i = 1;
      low_ip1 = 2;
      while (high_i > low_ip1) {
        mid_i = (low_i >> 1) + (high_i >> 1);
        if (((static_cast<unsigned int>(low_i) & 1U) == 1U) &&
            ((static_cast<unsigned int>(high_i) & 1U) == 1U)) {
          mid_i++;
        }
        if (x[ix] >= pp.breaks[mid_i - 1]) {
          low_i = mid_i;
          low_ip1 = mid_i + 1;
        } else {
          high_i = mid_i;
        }
      }
      xloc = x[ix] - pp.breaks[low_i - 1];
      b_v = pp.coefs[low_i - 1];
      for (high_i = 2; high_i <= numTerms; high_i++) {
        b_v = xloc * b_v + pp.coefs[(low_i + (high_i - 1) * coefStride) - 1];
      }
    }
    v[ix] = b_v;
  }
}

} // namespace coder

// End of code generation (ppval.cpp)
