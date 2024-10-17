//
// Trial License - for use to evaluate programs for possible purchase as
// an end-user only.
//
// besseli.cpp
//
// Code generation for function 'besseli'
//

// Include files
#include "besseli.h"
#include "casyi.h"
#include "cmlri.h"
#include "gammaln.h"
#include "log.h"
#include "rt_nonfinite.h"
#include <cmath>

// Function Definitions
namespace coder {
creal_T besseli(double z)
{
  creal_T hz;
  creal_T w;
  creal_T zd;
  double AZ;
  int inw;
  zd.re = z;
  zd.im = 0.0;
  if (std::isnan(z)) {
    w.re = rtNaN;
    w.im = 0.0;
  } else {
    double az;
    double b_atol;
    double im;
    double re;
    int ierr;
    int nw;
    boolean_T guard1;
    ierr = 0;
    b_atol = std::abs(z);
    if (b_atol > 0.0) {
      AZ = b_atol;
    } else {
      AZ = 0.0;
    }
    if (AZ > 1.0737418235E+9) {
      ierr = 4;
    } else if (AZ > 32767.999992370605) {
      ierr = 3;
    }
    w.re = 0.0;
    w.im = 0.0;
    if (b_atol > 0.0) {
      az = b_atol;
    } else {
      az = 0.0;
    }
    guard1 = false;
    if (az <= 2.0) {
      nw = 0;
      if (b_atol > 0.0) {
        double crsc_re;
        boolean_T iflag;
        crsc_re = 1.0;
        iflag = false;
        if (b_atol < 2.2250738585072014E-305) {
          w.re = 1.0;
          w.im = 0.0;
        } else {
          double acz;
          double cz_im;
          double cz_re;
          double s1_re;
          hz.re = 0.5 * z;
          hz.im = 0.0;
          if (b_atol > 4.7170688552396617E-153) {
            cz_re = hz.re * hz.re;
            AZ = hz.re * 0.0;
            cz_im = AZ + AZ;
            acz = std::abs(cz_re);
            if (!(acz > cz_im)) {
              if (std::isnan(cz_im)) {
                acz = rtNaN;
              } else {
                acz *= 1.4142135623730951;
              }
            }
          } else {
            cz_re = 0.0;
            cz_im = 0.0;
            acz = 0.0;
          }
          b_log(hz);
          AZ = 1.0;
          gammaln(AZ);
          s1_re = hz.re * 0.0 - AZ;
          AZ = hz.im * 0.0;
          if (s1_re > -700.92179369444591) {
            double aa;
            double ascle;
            double coef_re;
            double s1_im;
            boolean_T guard2;
            ascle = 0.0;
            if (s1_re <= -664.87164553371019) {
              iflag = true;
              crsc_re = 2.2204460492503131E-16;
              ascle = 1.0020841800044864E-289;
            }
            aa = std::exp(s1_re);
            if (iflag) {
              aa /= 2.2204460492503131E-16;
            }
            coef_re = aa * std::cos(AZ);
            AZ = aa * std::sin(AZ);
            b_atol = 2.2204460492503131E-16 * acz;
            s1_re = 1.0;
            s1_im = 0.0;
            if (!(acz < 2.2204460492503131E-16)) {
              double ak;
              double s;
              hz.re = 1.0;
              hz.im = 0.0;
              ak = 3.0;
              s = 1.0;
              aa = 2.0;
              double rs;
              do {
                rs = 1.0 / s;
                re = hz.re * cz_re - hz.im * cz_im;
                im = hz.re * cz_im + hz.im * cz_re;
                hz.re = rs * re;
                hz.im = rs * im;
                s1_re += hz.re;
                s1_im += hz.im;
                s += ak;
                ak += 2.0;
                aa = aa * acz * rs;
              } while (!!(aa > b_atol));
            }
            hz.re = s1_re * coef_re - s1_im * AZ;
            s1_re = s1_re * AZ + s1_im * coef_re;
            guard2 = false;
            if (iflag) {
              AZ = std::abs(hz.re);
              if (AZ > s1_re) {
                b_atol = 0.0;
              } else {
                b_atol = AZ;
                AZ = s1_re;
              }
              if ((!(b_atol <= ascle)) ||
                  (!(AZ < b_atol / 2.2204460492503131E-16))) {
                guard2 = true;
              }
            } else {
              guard2 = true;
            }
            if (guard2) {
              w.re = hz.re * crsc_re - s1_re * 0.0;
              w.im = hz.re * 0.0 + s1_re * crsc_re;
            }
          } else {
            nw = 1;
            if (acz > 0.0) {
              nw = -1;
            }
          }
        }
      } else {
        w.re = 1.0;
        w.im = 0.0;
      }
      if (nw < 0) {
        inw = 1;
      } else {
        inw = nw;
      }
      if ((1 - inw != 0) && (nw < 0)) {
        guard1 = true;
      }
    } else {
      guard1 = true;
    }
    if (guard1) {
      if (az < 21.784271729432426) {
        nw = cmlri(zd, w);
        if (nw < 0) {
          if (nw == -2) {
            inw = -2;
          } else {
            inw = -1;
          }
        } else {
          inw = 0;
        }
      } else {
        nw = casyi(zd, w);
        if (nw < 0) {
          if (nw == -2) {
            inw = -2;
          } else {
            inw = -1;
          }
        } else {
          inw = 0;
        }
      }
    }
    if (inw < 0) {
      if (inw == -2) {
        ierr = 5;
      } else {
        ierr = 2;
      }
    } else if ((!(z >= 0.0)) && (inw != 1)) {
      if (std::fmax(std::abs(w.re), std::abs(w.im)) <=
          1.0020841800044864E-289) {
        w.re *= 4.503599627370496E+15;
        w.im *= 4.503599627370496E+15;
        AZ = 2.2204460492503131E-16;
      } else {
        AZ = 1.0;
      }
      re = w.re - w.im * 0.0;
      im = w.re * 0.0 + w.im;
      w.re = AZ * re;
      w.im = AZ * im;
    }
    if (ierr == 5) {
      w.re = rtNaN;
      w.im = 0.0;
    } else if (ierr == 2) {
      w.re = rtInf;
      w.im = 0.0;
    }
    if (z > 0.0) {
      AZ = w.re;
      w.re = AZ;
      w.im = 0.0;
    }
  }
  return w;
}

} // namespace coder

// End of code generation (besseli.cpp)
