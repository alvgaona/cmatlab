//
// Trial License - for use to evaluate programs for possible purchase as
// an end-user only.
//
// fftSqueeze_terminate.cpp
//
// Code generation for function 'fftSqueeze_terminate'
//

// Include files
#include "fftSqueeze_terminate.h"
#include "fftSqueeze_data.h"
#include "rt_nonfinite.h"
#include "omp.h"

// Function Definitions
void fftSqueeze_terminate()
{
  omp_destroy_nest_lock(&fftSqueeze_nestLockGlobal);
  isInitialized_fftSqueeze = false;
}

// End of code generation (fftSqueeze_terminate.cpp)
