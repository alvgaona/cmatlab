//
// Trial License - for use to evaluate programs for possible purchase as
// an end-user only.
//
// resizeim_terminate.cpp
//
// Code generation for function 'resizeim_terminate'
//

// Include files
#include "resizeim_terminate.h"
#include "resizeim_data.h"
#include "rt_nonfinite.h"
#include "omp.h"

// Function Definitions
void resizeim_terminate()
{
  omp_destroy_nest_lock(&resizeim_nestLockGlobal);
  isInitialized_resizeim = false;
}

// End of code generation (resizeim_terminate.cpp)
