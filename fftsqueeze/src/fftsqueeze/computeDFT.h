//
// Trial License - for use to evaluate programs for possible purchase as
// an end-user only.
//
// computeDFT.h
//
// Code generation for function 'computeDFT'
//

#ifndef COMPUTEDFT_H
#define COMPUTEDFT_H

// Include files
#include "rtwtypes.h"
#include "coder_array.h"
#include <cstddef>
#include <cstdlib>

// Function Declarations
namespace coder {
void computeDFTviaFFT(const array<double, 2U> &xin, double nx, double nfft,
                      double Fs, array<creal_T, 2U> &Xx, array<double, 1U> &f);

}

#endif
// End of code generation (computeDFT.h)
