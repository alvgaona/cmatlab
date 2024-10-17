//
// Trial License - for use to evaluate programs for possible purchase as
// an end-user only.
//
// fftsqueeze.h
//
// Code generation for function 'fftsqueeze'
//

#ifndef FFTSQUEEZE_H
#define FFTSQUEEZE_H

// Include files
#include "rtwtypes.h"
#include "coder_array.h"
#include <cstddef>
#include <cstdlib>

// Function Declarations
extern void fftsqueeze(const coder::array<double, 2U> &x, double fs,
                       const coder::array<double, 1U> &window,
                       coder::array<creal_T, 2U> &s,
                       coder::array<double, 1U> &f,
                       coder::array<double, 2U> &t);

#endif
// End of code generation (fftsqueeze.h)
