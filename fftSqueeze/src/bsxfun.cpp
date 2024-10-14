//
// Trial License - for use to evaluate programs for possible purchase as
// an end-user only.
//
// bsxfun.cpp
//
// Code generation for function 'bsxfun'
//

// Include files
#include "bsxfun.h"
#include "rt_nonfinite.h"
#include "coder_array.h"

// Function Definitions
namespace coder {
void bsxfun(const double a_data[], int a_size, const array<double, 2U> &b,
            array<double, 2U> &c)
{
  int i;
  int i1;
  if (a_size == 1) {
    i = 128;
  } else if (a_size == 128) {
    i = 128;
  } else {
    i = a_size;
  }
  i1 = b.size(1);
  c.set_size(i, b.size(1));
  if (b.size(1) != 0) {
    int acoef;
    int bcoef;
    bcoef = (b.size(1) != 1);
    acoef = (a_size != 1);
    for (int k{0}; k < i1; k++) {
      int varargin_3;
      varargin_3 = bcoef * k;
      for (int b_k{0}; b_k < i; b_k++) {
        c[b_k + c.size(0) * k] =
            a_data[acoef * b_k] * b[b_k + 128 * varargin_3];
      }
    }
  }
}

} // namespace coder

// End of code generation (bsxfun.cpp)
