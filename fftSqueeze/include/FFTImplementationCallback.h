//
// Trial License - for use to evaluate programs for possible purchase as
// an end-user only.
//
// FFTImplementationCallback.h
//
// Code generation for function 'FFTImplementationCallback'
//

#ifndef FFTIMPLEMENTATIONCALLBACK_H
#define FFTIMPLEMENTATIONCALLBACK_H

// Include files
#include "rtwtypes.h"
#include "coder_array.h"
#include <cstddef>
#include <cstdlib>

// Type Definitions
namespace coder {
namespace internal {
namespace fft {
class FFTImplementationCallback {
public:
  static void r2br_r2dit_trig(const array<double, 2U> &x,
                              const double costab_data[],
                              const double sintab_data[],
                              array<creal_T, 2U> &y);

protected:
  static void doHalfLengthRadix2(const array<double, 2U> &x, int xoffInit,
                                 creal_T y_data[], const double costab_data[],
                                 const double sintab_data[]);
};

} // namespace fft
} // namespace internal
} // namespace coder

#endif
// End of code generation (FFTImplementationCallback.h)
