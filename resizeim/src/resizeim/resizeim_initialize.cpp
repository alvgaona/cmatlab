//
// Trial License - for use to evaluate programs for possible purchase as
// an end-user only.
//
// resizeim_initialize.cpp
//
// Code generation for function 'resizeim_initialize'
//

// Include files
#include "resizeim_initialize.h"
#include "resizeim_data.h"
#include "rt_nonfinite.h"
#include "omp.h"

// Function Definitions
void resizeim_initialize()
{
  omp_init_nest_lock(&resizeim_nestLockGlobal);
  isInitialized_resizeim = true;
}

// End of code generation (resizeim_initialize.cpp)
