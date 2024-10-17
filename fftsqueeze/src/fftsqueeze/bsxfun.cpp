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
void bsxfun(const array<double, 1U> &a, const array<double, 2U> &b,
            array<double, 2U> &c)
{
  int u0;
  int u1;
  u0 = b.size(0);
  u1 = a.size(0);
  if (u0 <= u1) {
    u1 = u0;
  }
  if (b.size(0) == 1) {
    u0 = a.size(0);
  } else if (a.size(0) == 1) {
    u0 = b.size(0);
  } else if (a.size(0) == b.size(0)) {
    u0 = a.size(0);
  } else {
    u0 = u1;
  }
  u1 = b.size(1);
  c.set_size(u0, b.size(1));
  if ((u0 != 0) && (b.size(1) != 0)) {
    int acoef;
    int b_bcoef;
    int bcoef;
    bcoef = (b.size(1) != 1);
    acoef = (a.size(0) != 1);
    b_bcoef = (b.size(0) != 1);
    for (int k{0}; k < u1; k++) {
      int varargin_3;
      varargin_3 = bcoef * k;
      for (int b_k{0}; b_k < u0; b_k++) {
        c[b_k + c.size(0) * k] =
            a[acoef * b_k] * b[b_bcoef * b_k + b.size(0) * varargin_3];
      }
    }
  }
}

} // namespace coder

// End of code generation (bsxfun.cpp)
