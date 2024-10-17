//
// Trial License - for use to evaluate programs for possible purchase as
// an end-user only.
//
// fftsqueeze_initialize.cpp
//
// Code generation for function 'fftsqueeze_initialize'
//

// Include files
#include "fftsqueeze_initialize.h"
#include "fftsqueeze_data.h"
#include "rt_nonfinite.h"
#include "omp.h"

// Function Definitions
void fftsqueeze_initialize()
{
  omp_init_nest_lock(&fftsqueeze_nestLockGlobal);
  isInitialized_fftsqueeze = true;
}

// End of code generation (fftsqueeze_initialize.cpp)
