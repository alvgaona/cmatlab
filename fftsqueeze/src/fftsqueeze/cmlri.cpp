//
// Trial License - for use to evaluate programs for possible purchase as
// an end-user only.
//
// cmlri.cpp
//
// Code generation for function 'cmlri'
//

// Include files
#include "cmlri.h"
#include "fftsqueeze_data.h"
#include "gammaln.h"
#include "log.h"
#include "rt_nonfinite.h"
#include <cmath>

// Function Definitions
namespace coder {
int cmlri(const creal_T z, creal_T &y)
{
  creal_T rz;
  double ack;
  double ak;
  double ap;
  double az;
  double b_tmp;
  double ck_im;
  double ck_re;
  double fkk;
  double flooraz;
  double p1_im;
  double p1_re;
  double p2_im;
  double p2_re;
  double pt_im;
  double pt_re;
  double rho;
  double rho2;
  double tst;
  int i;
  int icounter;
  int nz;
  boolean_T earlyExit;
  boolean_T exitg1;
  nz = 0;
  fkk = std::abs(z.re);
  b_tmp = std::abs(z.im);
  if (fkk < b_tmp) {
    rho = fkk / b_tmp;
    az = b_tmp * std::sqrt(rho * rho + 1.0);
  } else if (fkk > b_tmp) {
    rho2 = b_tmp / fkk;
    az = fkk * std::sqrt(rho2 * rho2 + 1.0);
  } else if (std::isnan(b_tmp)) {
    az = rtNaN;
  } else {
    az = fkk * 1.4142135623730951;
  }
  flooraz = std::floor(az);
  if (z.im == 0.0) {
    ck_re = (flooraz + 1.0) / z.re;
    ck_im = 0.0;
    rz.re = 2.0 / z.re;
    rz.im = 0.0;
  } else if (z.re == 0.0) {
    ck_re = 0.0;
    ck_im = -((flooraz + 1.0) / z.im);
    rz.re = 0.0;
    rz.im = -(2.0 / z.im);
  } else if (fkk > b_tmp) {
    rho2 = z.im / z.re;
    rho = z.re + rho2 * z.im;
    ck_re = ((flooraz + 1.0) + rho2 * 0.0) / rho;
    ck_im = (0.0 - rho2 * (flooraz + 1.0)) / rho;
    rho = z.re + rho2 * z.im;
    rz.re = (rho2 * 0.0 + 2.0) / rho;
    rz.im = (0.0 - rho2 * 2.0) / rho;
  } else if (b_tmp == fkk) {
    if (z.re > 0.0) {
      rho2 = 0.5;
    } else {
      rho2 = -0.5;
    }
    if (z.im > 0.0) {
      rho = 0.5;
    } else {
      rho = -0.5;
    }
    ck_re = ((flooraz + 1.0) * rho2 + 0.0 * rho) / fkk;
    ck_im = (0.0 * rho2 - (flooraz + 1.0) * rho) / fkk;
    if (z.re > 0.0) {
      rho2 = 0.5;
    } else {
      rho2 = -0.5;
    }
    if (z.im > 0.0) {
      rho = 0.5;
    } else {
      rho = -0.5;
    }
    rz.re = (2.0 * rho2 + 0.0 * rho) / fkk;
    rz.im = (0.0 * rho2 - 2.0 * rho) / fkk;
  } else {
    rho2 = z.re / z.im;
    rho = z.im + rho2 * z.re;
    ck_re = rho2 * (flooraz + 1.0) / rho;
    ck_im = (rho2 * 0.0 - (flooraz + 1.0)) / rho;
    rho = z.im + rho2 * z.re;
    rz.re = rho2 * 2.0 / rho;
    rz.im = (rho2 * 0.0 - 2.0) / rho;
  }
  p1_re = 0.0;
  p1_im = 0.0;
  p2_re = 1.0;
  p2_im = 0.0;
  ack = ((flooraz + 1.0) + 1.0) / az;
  rho = ack + std::sqrt(ack * ack - 1.0);
  rho2 = rho * rho;
  tst = (rho2 + rho2) / ((rho2 - 1.0) * (rho - 1.0)) / 2.2204460492503131E-16;
  ak = flooraz + 1.0;
  earlyExit = true;
  icounter = 1;
  i = 0;
  exitg1 = false;
  while ((!exitg1) && (i < 80)) {
    icounter++;
    pt_re = p2_re;
    pt_im = p2_im;
    rho = ck_re * p2_re - ck_im * p2_im;
    rho2 = ck_re * p2_im + ck_im * p2_re;
    p2_re = p1_re - rho;
    p2_im = p1_im - rho2;
    p1_re = pt_re;
    p1_im = pt_im;
    ck_re += rz.re;
    ck_im += rz.im;
    rho = std::abs(p2_re);
    rho2 = std::abs(p2_im);
    if (rho < rho2) {
      rho /= rho2;
      ap = rho2 * std::sqrt(rho * rho + 1.0);
    } else if (rho > rho2) {
      rho2 /= rho;
      ap = rho * std::sqrt(rho2 * rho2 + 1.0);
    } else if (std::isnan(rho2)) {
      ap = rtNaN;
    } else {
      ap = rho * 1.4142135623730951;
    }
    if (ap > tst * ak * ak) {
      earlyExit = false;
      exitg1 = true;
    } else {
      ak++;
      i++;
    }
  }
  if (earlyExit) {
    nz = -2;
  } else {
    int kcounter;
    boolean_T guard1;
    kcounter = 1;
    guard1 = false;
    if (static_cast<int>(flooraz) <= 0) {
      int itime;
      p1_re = 0.0;
      p1_im = 0.0;
      p2_re = 1.0;
      p2_im = 0.0;
      if (z.im == 0.0) {
        ck_re = 1.0 / z.re;
        ck_im = 0.0;
      } else if (z.re == 0.0) {
        ck_re = 0.0;
        ck_im = -(1.0 / z.im);
      } else if (fkk > b_tmp) {
        rho2 = z.im / z.re;
        rho = z.re + rho2 * z.im;
        ck_re = (rho2 * 0.0 + 1.0) / rho;
        ck_im = (0.0 - rho2) / rho;
      } else if (b_tmp == fkk) {
        if (z.re > 0.0) {
          rho2 = 0.5;
        } else {
          rho2 = -0.5;
        }
        if (z.im > 0.0) {
          rho = 0.5;
        } else {
          rho = -0.5;
        }
        ck_re = (rho2 + 0.0 * rho) / fkk;
        ck_im = (0.0 * rho2 - rho) / fkk;
      } else {
        rho2 = z.re / z.im;
        rho = z.im + rho2 * z.re;
        ck_re = rho2 / rho;
        ck_im = (rho2 * 0.0 - 1.0) / rho;
      }
      tst = std::sqrt(1.0 / az / 2.2204460492503131E-16);
      itime = 1;
      earlyExit = true;
      i = 0;
      exitg1 = false;
      while ((!exitg1) && (i < 80)) {
        kcounter++;
        pt_re = p2_re;
        pt_im = p2_im;
        rho = ck_re * p2_re - ck_im * p2_im;
        rho2 = ck_re * p2_im + ck_im * p2_re;
        p2_re = p1_re - rho;
        p2_im = p1_im - rho2;
        p1_re = pt_re;
        p1_im = pt_im;
        ck_re += rz.re;
        ck_im += rz.im;
        rho = std::abs(p2_re);
        rho2 = std::abs(p2_im);
        if (rho < rho2) {
          rho /= rho2;
          ap = rho2 * std::sqrt(rho * rho + 1.0);
        } else if (rho > rho2) {
          rho2 /= rho;
          ap = rho * std::sqrt(rho2 * rho2 + 1.0);
        } else if (std::isnan(rho2)) {
          ap = rtNaN;
        } else {
          ap = rho * 1.4142135623730951;
        }
        if (ap >= tst * ak * ak) {
          if (itime == 2) {
            earlyExit = false;
            exitg1 = true;
          } else {
            rho = std::abs(ck_re);
            rho2 = std::abs(ck_im);
            if (rho < rho2) {
              rho /= rho2;
              ack = rho2 * std::sqrt(rho * rho + 1.0);
            } else if (rho > rho2) {
              rho2 /= rho;
              ack = rho * std::sqrt(rho2 * rho2 + 1.0);
            } else if (std::isnan(rho2)) {
              ack = rtNaN;
            } else {
              ack = rho * 1.4142135623730951;
            }
            rho = std::abs(pt_re);
            rho2 = std::abs(pt_im);
            if (rho < rho2) {
              rho /= rho2;
              rho = rho2 * std::sqrt(rho * rho + 1.0);
            } else if (rho > rho2) {
              rho2 /= rho;
              rho *= std::sqrt(rho2 * rho2 + 1.0);
            } else if (std::isnan(rho2)) {
              rho = rtNaN;
            } else {
              rho *= 1.4142135623730951;
            }
            rho = std::fmin(ack + std::sqrt(ack * ack - 1.0), ap / rho);
            tst *= std::sqrt(rho / (rho * rho - 1.0));
            itime = 2;
            i++;
          }
        } else {
          i++;
        }
      }
      if (earlyExit) {
        nz = -2;
      } else {
        guard1 = true;
      }
    } else {
      guard1 = true;
    }
    if (guard1) {
      icounter += static_cast<int>(flooraz);
      if (icounter >= kcounter) {
        kcounter = icounter;
      }
      fkk = kcounter;
      p1_re = 0.0;
      p1_im = 0.0;
      p2_re = 1.0020841800044864E-289;
      p2_im = 0.0;
      rho2 = static_cast<double>(kcounter) + 1.0;
      gammaln(rho2);
      b_tmp = 1.0;
      gammaln(b_tmp);
      tst = std::exp((rho2 - rho2) - b_tmp);
      ak = 0.0;
      az = 0.0;
      for (i = 0; i < kcounter; i++) {
        pt_re = p2_re;
        pt_im = p2_im;
        rho2 = fkk * rz.re;
        rho = fkk * rz.im;
        ap = rho2 * p2_re - rho * p2_im;
        rho = rho2 * p2_im + rho * p2_re;
        p2_re = p1_re + ap;
        p2_im = p1_im + rho;
        p1_re = pt_re;
        p1_im = pt_im;
        ack = tst * (1.0 - 0.0 / fkk);
        rho2 = ack + tst;
        ak += rho2 * pt_re;
        az += rho2 * pt_im;
        tst = ack;
        fkk--;
      }
      y.re = p2_re;
      y.im = p2_im;
      b_log(rz);
      rho2 = 0.0 * rz.im;
      rho = 0.0 * rz.re;
      ck_re = ((rho - rho2) + z.re) - b_tmp;
      ck_im = (rho2 + rho) + z.im;
      p2_re += ak;
      p2_im += az;
      rho = std::abs(p2_re);
      rho2 = std::abs(p2_im);
      if (rho < rho2) {
        rho /= rho2;
        ap = rho2 * std::sqrt(rho * rho + 1.0);
      } else if (rho > rho2) {
        rho2 /= rho;
        ap = rho * std::sqrt(rho2 * rho2 + 1.0);
      } else if (std::isnan(rho2)) {
        ap = rtNaN;
      } else {
        ap = rho * 1.4142135623730951;
      }
      p1_re = 1.0 / ap;
      if (ck_re == 0.0) {
        ck_re = std::cos(ck_im);
        ck_im = std::sin(ck_im);
      } else if (ck_im == 0.0) {
        ck_re = std::exp(ck_re);
        ck_im = 0.0;
      } else if (std::isinf(ck_im) && std::isinf(ck_re) && (ck_re < 0.0)) {
        ck_re = 0.0;
        ck_im = 0.0;
      } else {
        rho2 = std::exp(ck_re / 2.0);
        ck_re = rho2 * (rho2 * std::cos(ck_im));
        ck_im = rho2 * (rho2 * std::sin(ck_im));
      }
      rho = ck_re * p1_re - ck_im * 0.0;
      ck_im = ck_re * 0.0 + ck_im * p1_re;
      rho2 = p2_re * p1_re + p2_im * 0.0;
      p2_im = p2_re * 0.0 - p2_im * p1_re;
      ck_re = rho * rho2 - ck_im * p2_im;
      ck_im = rho * p2_im + ck_im * rho2;
      rho2 = y.re * ck_im + y.im * ck_re;
      y.re = y.re * ck_re - y.im * ck_im;
      y.im = rho2;
    }
  }
  return nz;
}

} // namespace coder

// End of code generation (cmlri.cpp)
