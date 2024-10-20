//
// Trial License - for use to evaluate programs for possible purchase as
// an end-user only.
//
// rt_nonfinite.cpp
//
// Code generation for function 'fftsqueeze'
//

// Abstract:
//      MATLAB for code generation function to initialize non-finites,
//      (Inf, NaN and -Inf).
// Include files
#include "rt_nonfinite.h"
#include <cmath>
#include <limits>

real_T rtNaN{std::numeric_limits<real_T>::quiet_NaN()};
real_T rtInf{std::numeric_limits<real_T>::infinity()};
real_T rtMinusInf{-std::numeric_limits<real_T>::infinity()};
real32_T rtNaNF{std::numeric_limits<real32_T>::quiet_NaN()};
real32_T rtInfF{std::numeric_limits<real32_T>::infinity()};
real32_T rtMinusInfF{-std::numeric_limits<real32_T>::infinity()};

// End of code generation (rt_nonfinite.cpp)
