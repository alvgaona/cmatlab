//
// Trial License - for use to evaluate programs for possible purchase as
// an end-user only.
//
// ppval.h
//
// Code generation for function 'ppval'
//

#ifndef PPVAL_H
#define PPVAL_H

// Include files
#include "rtwtypes.h"
#include "coder_array.h"
#include <cstddef>
#include <cstdlib>

// Type Declarations
struct struct_T;

// Function Declarations
namespace coder {
void ppval(const struct_T &pp, const array<double, 1U> &x,
           array<double, 1U> &v);

}

#endif
// End of code generation (ppval.h)
