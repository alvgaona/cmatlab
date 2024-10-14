//
// Trial License - for use to evaluate programs for possible purchase as
// an end-user only.
//
// fftSqueeze_internal_types.h
//
// Code generation for function 'fftSqueeze'
//

#ifndef FFTSQUEEZE_INTERNAL_TYPES_H
#define FFTSQUEEZE_INTERNAL_TYPES_H

// Include files
#include "fftSqueeze_types.h"
#include "rtwtypes.h"
#include "coder_bounded_array.h"

// Type Definitions
struct struct_T {
  coder::bounded_array<double, 128U, 2U> breaks;
  coder::bounded_array<double, 508U, 2U> coefs;
};

#endif
// End of code generation (fftSqueeze_internal_types.h)
