//
// Trial License - for use to evaluate programs for possible purchase as
// an end-user only.
//
// casyi.cpp
//
// Code generation for function 'casyi'
//

// Include files
#include "casyi.h"
#include "fftSqueeze_data.h"
#include "rt_nonfinite.h"
#include <cmath>

// Function Definitions
namespace coder {
int casyi(const creal_T z, creal_T &y)
{
  double absxi;
  double absxr;
  double ak1_im;
  double ak1_re;
  double az;
  double r;
  double sgn;
  double yr;
  int nz;
  nz = 0;
  sgn = std::abs(z.re);
  absxi = std::abs(z.im);
  if (sgn < absxi) {
    r = sgn / absxi;
    az = absxi * std::sqrt(r * r + 1.0);
  } else if (sgn > absxi) {
    yr = absxi / sgn;
    az = sgn * std::sqrt(yr * yr + 1.0);
  } else if (std::isnan(absxi)) {
    az = rtNaN;
  } else {
    az = sgn * 1.4142135623730951;
  }
  if (z.im == 0.0) {
    ak1_re = 0.15915494309189535 / z.re;
    ak1_im = 0.0;
  } else if (z.re == 0.0) {
    ak1_re = 0.0;
    ak1_im = -(0.15915494309189535 / z.im);
  } else if (sgn > absxi) {
    absxi = z.im / z.re;
    absxr = z.re + absxi * z.im;
    ak1_re = (absxi * 0.0 + 0.15915494309189535) / absxr;
    ak1_im = (0.0 - absxi * 0.15915494309189535) / absxr;
  } else if (absxi == sgn) {
    if (z.re > 0.0) {
      absxi = 0.5;
    } else {
      absxi = -0.5;
    }
    if (z.im > 0.0) {
      yr = 0.5;
    } else {
      yr = -0.5;
    }
    ak1_re = (0.15915494309189535 * absxi + 0.0 * yr) / sgn;
    ak1_im = (0.0 * absxi - 0.15915494309189535 * yr) / sgn;
  } else {
    absxi = z.re / z.im;
    absxr = z.im + absxi * z.re;
    ak1_re = absxi * 0.15915494309189535 / absxr;
    ak1_im = (absxi * 0.0 - 0.15915494309189535) / absxr;
  }
  if (ak1_im == 0.0) {
    if (ak1_re < 0.0) {
      yr = 0.0;
      absxi = std::sqrt(-ak1_re);
    } else {
      yr = std::sqrt(ak1_re);
      absxi = 0.0;
    }
  } else if (ak1_re == 0.0) {
    if (ak1_im < 0.0) {
      yr = std::sqrt(-ak1_im / 2.0);
      absxi = -yr;
    } else {
      yr = std::sqrt(ak1_im / 2.0);
      absxi = yr;
    }
  } else if (std::isnan(ak1_re)) {
    yr = rtNaN;
    absxi = rtNaN;
  } else if (std::isnan(ak1_im)) {
    yr = rtNaN;
    absxi = rtNaN;
  } else if (std::isinf(ak1_im)) {
    yr = std::abs(ak1_im);
    absxi = ak1_im;
  } else if (std::isinf(ak1_re)) {
    if (ak1_re < 0.0) {
      yr = 0.0;
      absxi = ak1_im * -ak1_re;
    } else {
      yr = ak1_re;
      absxi = 0.0;
    }
  } else {
    absxr = std::abs(ak1_re);
    absxi = std::abs(ak1_im);
    if ((absxr > 4.4942328371557893E+307) ||
        (absxi > 4.4942328371557893E+307)) {
      absxr *= 0.5;
      yr = absxi * 0.5;
      if (absxr < yr) {
        r = absxr / yr;
        absxi = yr * std::sqrt(r * r + 1.0);
      } else if (absxr > yr) {
        yr /= absxr;
        absxi = absxr * std::sqrt(yr * yr + 1.0);
      } else {
        absxi = absxr * 1.4142135623730951;
      }
      if (absxi > absxr) {
        yr = std::sqrt(absxi) * std::sqrt(absxr / absxi + 1.0);
      } else {
        yr = std::sqrt(absxi) * 1.4142135623730951;
      }
    } else {
      if (absxr < absxi) {
        r = absxr / absxi;
        absxi *= std::sqrt(r * r + 1.0);
      } else if (absxr > absxi) {
        yr = absxi / absxr;
        absxi = absxr * std::sqrt(yr * yr + 1.0);
      } else {
        absxi = absxr * 1.4142135623730951;
      }
      yr = std::sqrt((absxi + absxr) * 0.5);
    }
    if (ak1_re > 0.0) {
      absxi = 0.5 * (ak1_im / yr);
    } else {
      if (ak1_im < 0.0) {
        absxi = -yr;
      } else {
        absxi = yr;
      }
      yr = 0.5 * (ak1_im / absxi);
    }
  }
  if (sgn > 700.92179369444591) {
    nz = -1;
    y.re = rtNaN;
    y.im = 0.0;
  } else {
    double aa;
    double b_re;
    double bb;
    double cs1_im;
    double cs1_re;
    double cs2_im;
    double cs2_re;
    double dk_im;
    double dk_re;
    double ez_im;
    double ez_re;
    double im;
    double re;
    double tmp_im;
    double tmp_re;
    int bk;
    signed char p1_im;
    boolean_T errflag;
    boolean_T exitg1;
    if (z.re == 0.0) {
      tmp_re = std::cos(z.im);
      tmp_im = std::sin(z.im);
    } else if (z.im == 0.0) {
      tmp_re = std::exp(z.re);
      tmp_im = 0.0;
    } else if (std::isinf(z.im) && std::isinf(z.re) && (z.re < 0.0)) {
      tmp_re = 0.0;
      tmp_im = 0.0;
    } else {
      r = std::exp(z.re / 2.0);
      tmp_re = r * (r * std::cos(z.im));
      tmp_im = r * (r * std::sin(z.im));
    }
    re = yr * tmp_re - absxi * tmp_im;
    im = yr * tmp_im + absxi * tmp_re;
    ez_re = 8.0 * z.re;
    ez_im = 8.0 * z.im;
    ak1_im = 8.0 * az;
    if (z.im != 0.0) {
      bk = 1;
      if (z.im < 0.0) {
        bk = -1;
      }
      p1_im = static_cast<signed char>(bk);
    } else {
      p1_im = 0;
    }
    r = -1.0;
    ak1_re = 2.2204460492503131E-16 / ak1_im;
    sgn = 1.0;
    cs1_re = 1.0;
    cs1_im = 0.0;
    cs2_re = 1.0;
    cs2_im = 0.0;
    tmp_re = 1.0;
    tmp_im = 0.0;
    az = 0.0;
    aa = 1.0;
    bb = ak1_im;
    dk_re = ez_re;
    dk_im = ez_im;
    errflag = true;
    bk = 0;
    exitg1 = false;
    while ((!exitg1) && (bk < 45)) {
      tmp_re *= r;
      tmp_im *= r;
      if (dk_im == 0.0) {
        if (tmp_im == 0.0) {
          b_re = tmp_re / dk_re;
          tmp_im = 0.0;
        } else if (tmp_re == 0.0) {
          b_re = 0.0;
          tmp_im /= dk_re;
        } else {
          b_re = tmp_re / dk_re;
          tmp_im /= dk_re;
        }
      } else if (dk_re == 0.0) {
        if (tmp_re == 0.0) {
          b_re = tmp_im / dk_im;
          tmp_im = 0.0;
        } else if (tmp_im == 0.0) {
          b_re = 0.0;
          tmp_im = -(tmp_re / dk_im);
        } else {
          b_re = tmp_im / dk_im;
          tmp_im = -(tmp_re / dk_im);
        }
      } else {
        absxr = std::abs(dk_re);
        absxi = std::abs(dk_im);
        if (absxr > absxi) {
          absxi = dk_im / dk_re;
          absxr = dk_re + absxi * dk_im;
          b_re = (tmp_re + absxi * tmp_im) / absxr;
          tmp_im = (tmp_im - absxi * tmp_re) / absxr;
        } else if (absxi == absxr) {
          if (dk_re > 0.0) {
            absxi = 0.5;
          } else {
            absxi = -0.5;
          }
          if (dk_im > 0.0) {
            yr = 0.5;
          } else {
            yr = -0.5;
          }
          b_re = (tmp_re * absxi + tmp_im * yr) / absxr;
          tmp_im = (tmp_im * absxi - tmp_re * yr) / absxr;
        } else {
          absxi = dk_re / dk_im;
          absxr = dk_im + absxi * dk_re;
          b_re = (absxi * tmp_re + tmp_im) / absxr;
          tmp_im = (absxi * tmp_im - tmp_re) / absxr;
        }
      }
      tmp_re = b_re;
      cs2_re += b_re;
      cs2_im += tmp_im;
      sgn = -sgn;
      cs1_re += b_re * sgn;
      cs1_im += tmp_im * sgn;
      dk_re += ez_re;
      dk_im += ez_im;
      aa = aa * std::abs(r) / bb;
      bb += ak1_im;
      az += 8.0;
      r -= az;
      if (aa <= ak1_re) {
        errflag = false;
        exitg1 = true;
      } else {
        bk++;
      }
    }
    if (errflag) {
      nz = -2;
    } else {
      if (z.re + z.re < 700.92179369444591) {
        tmp_re = -2.0 * z.re;
        tmp_im = -2.0 * z.im;
        if (tmp_re == 0.0) {
          tmp_re = std::cos(tmp_im);
          tmp_im = std::sin(tmp_im);
        } else if (tmp_im == 0.0) {
          tmp_re = std::exp(tmp_re);
          tmp_im = 0.0;
        } else if (std::isinf(tmp_im) && std::isinf(tmp_re) && (tmp_re < 0.0)) {
          tmp_re = 0.0;
          tmp_im = 0.0;
        } else {
          r = std::exp(tmp_re / 2.0);
          tmp_re = r * (r * std::cos(tmp_im));
          tmp_im = r * (r * std::sin(tmp_im));
        }
        b_re = tmp_re * cs2_re - tmp_im * cs2_im;
        absxi = tmp_re * cs2_im + tmp_im * cs2_re;
        cs1_re += b_re * 0.0 - absxi * static_cast<double>(p1_im);
        cs1_im += b_re * static_cast<double>(p1_im) + absxi * 0.0;
      }
      y.re = cs1_re * re - cs1_im * im;
      y.im = cs1_re * im + cs1_im * re;
    }
  }
  return nz;
}

} // namespace coder

// End of code generation (casyi.cpp)
