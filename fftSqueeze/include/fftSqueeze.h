//
// Trial License - for use to evaluate programs for possible purchase as
// an end-user only.
//
// fftSqueeze.h
//
// Code generation for function 'fftSqueeze'
//

#ifndef FFTSQUEEZE_H
#define FFTSQUEEZE_H

// Include files
#include "rtwtypes.h"
#include "coder_array.h"
#include <cstddef>
#include <cstdlib>

// Function Declarations
extern void fftSqueeze(const coder::array<double, 2U> &x, double fs,
                       const double window[128], coder::array<creal_T, 2U> &s,
                       double f_data[], int f_size[1],
                       coder::array<double, 2U> &t);

#endif
// End of code generation (fftSqueeze.h)
