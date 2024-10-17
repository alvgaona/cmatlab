//
// Trial License - for use to evaluate programs for possible purchase as
// an end-user only.
//
// dtwin.h
//
// Code generation for function 'dtwin'
//

#ifndef DTWIN_H
#define DTWIN_H

// Include files
#include "rtwtypes.h"
#include "coder_array.h"
#include <cstddef>
#include <cstdlib>

// Function Declarations
namespace coder {
namespace b_signal {
namespace internal {
namespace spectral {
void dtwin(const array<double, 1U> &w, double Fs, array<double, 1U> &Wdt);

}
} // namespace internal
} // namespace b_signal
} // namespace coder

#endif
// End of code generation (dtwin.h)
