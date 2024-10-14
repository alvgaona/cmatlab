//
// Trial License - for use to evaluate programs for possible purchase as
// an end-user only.
//
// fsstParser.h
//
// Code generation for function 'fsstParser'
//

#ifndef FSSTPARSER_H
#define FSSTPARSER_H

// Include files
#include "rtwtypes.h"
#include "coder_array.h"
#include <cstddef>
#include <cstdlib>

// Function Declarations
namespace coder {
namespace b_signal {
namespace internal {
namespace fsst {
double fsstParser(const array<double, 2U> &x, double varargin_1,
                  const double varargin_2[128], double win_data[],
                  int &win_size);

}
} // namespace internal
} // namespace b_signal
} // namespace coder

#endif
// End of code generation (fsstParser.h)
