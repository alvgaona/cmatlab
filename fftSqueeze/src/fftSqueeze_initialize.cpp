//
// Trial License - for use to evaluate programs for possible purchase as
// an end-user only.
//
// fftSqueeze_initialize.cpp
//
// Code generation for function 'fftSqueeze_initialize'
//

// Include files
#include "fftSqueeze_initialize.h"
#include "fftSqueeze_data.h"
#include "rt_nonfinite.h"
#include "omp.h"

// Function Definitions
void fftSqueeze_initialize()
{
  omp_init_nest_lock(&fftSqueeze_nestLockGlobal);
  isInitialized_fftSqueeze = true;
}

// End of code generation (fftSqueeze_initialize.cpp)
