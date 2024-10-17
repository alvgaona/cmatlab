//
// Trial License - for use to evaluate programs for possible purchase as
// an end-user only.
//
// fftsqueeze_terminate.cpp
//
// Code generation for function 'fftsqueeze_terminate'
//

// Include files
#include "fftsqueeze_terminate.h"
#include "fftsqueeze_data.h"
#include "rt_nonfinite.h"
#include "omp.h"

// Function Definitions
void fftsqueeze_terminate()
{
  omp_destroy_nest_lock(&fftsqueeze_nestLockGlobal);
  isInitialized_fftsqueeze = false;
}

// End of code generation (fftsqueeze_terminate.cpp)
