//
// Trial License - for use to evaluate programs for possible purchase as
// an end-user only.
//
// log.cpp
//
// Code generation for function 'log'
//

// Include files
#include "log.h"
#include "atan2.h"
#include "rt_nonfinite.h"
#include <cmath>

// Function Definitions
namespace coder {
void b_log(creal_T &x)
{
  if (x.im == 0.0) {
    if (x.re < 0.0) {
      x.re = std::log(std::abs(x.re));
      x.im = 3.1415926535897931;
    } else {
      x.re = std::log(std::abs(x.re));
      x.im = 0.0;
    }
  } else {
    double a;
    double b;
    boolean_T guard1;
    a = std::abs(x.re);
    guard1 = false;
    if (a > 8.9884656743115785E+307) {
      guard1 = true;
    } else {
      b = std::abs(x.im);
      if (b > 8.9884656743115785E+307) {
        guard1 = true;
      } else {
        if (a < b) {
          a /= b;
          a = b * std::sqrt(a * a + 1.0);
        } else if (a > b) {
          b /= a;
          a *= std::sqrt(b * b + 1.0);
        } else if (std::isnan(b)) {
          a = rtNaN;
        } else {
          a *= 1.4142135623730951;
        }
        b = x.re;
        x.re = std::log(a);
        x.im = internal::scalar::b_atan2(x.im, b);
      }
    }
    if (guard1) {
      a = std::abs(x.re / 2.0);
      b = std::abs(x.im / 2.0);
      if (a < b) {
        a /= b;
        a = b * std::sqrt(a * a + 1.0);
      } else if (a > b) {
        b /= a;
        a *= std::sqrt(b * b + 1.0);
      } else if (std::isnan(b)) {
        a = rtNaN;
      } else {
        a *= 1.4142135623730951;
      }
      b = x.re;
      x.re = std::log(a) + 0.69314718055994529;
      x.im = internal::scalar::b_atan2(x.im, b);
    }
  }
}

} // namespace coder

// End of code generation (log.cpp)
