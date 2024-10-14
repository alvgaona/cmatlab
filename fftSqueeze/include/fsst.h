//
// Trial License - for use to evaluate programs for possible purchase as
// an end-user only.
//
// fsst.h
//
// Code generation for function 'fsst'
//

#ifndef FSST_H
#define FSST_H

// Include files
#include "rtwtypes.h"
#include "coder_array.h"
#include <cstddef>
#include <cstdlib>

// Function Declarations
namespace coder {
int fsst(const array<double, 2U> &x, double varargin_1,
         const double varargin_2[128], array<creal_T, 2U> &sst, double f_data[],
         array<double, 2U> &t);

}

#endif
// End of code generation (fsst.h)
