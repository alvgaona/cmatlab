
//
// Trial License - for use to evaluate programs for possible purchase as
// an end-user only.
//
// main.cpp
//
// Code generation for function 'main'
//

/*************************************************************************/
/* This automatically generated example C++ main file shows how to call  */
/* entry-point functions that MATLAB Coder generated. You must customize */
/* this file for your application. Do not modify this file directly.     */
/* Instead, make a copy of this file, modify it, and integrate it into   */
/* your development environment.                                         */
/*                                                                       */
/* This file initializes entry-point function arguments to a default     */
/* size and value before calling the entry-point functions. It does      */
/* not store or use any values returned from the entry-point functions.  */
/* If necessary, it does pre-allocate memory for returned values.        */
/* You can use this file as a starting point for a main function that    */
/* you can deploy in your application.                                   */
/*                                                                       */
/* After you copy the file, and before you deploy it, you must make the  */
/* following changes:                                                    */
/* * For variable-size function arguments, change the example sizes to   */
/* the sizes that your application requires.                             */
/* * Change the example values of function arguments to the values that  */
/* your application requires.                                            */
/* * If the entry-point functions return values, store these values or   */
/* otherwise use them as required by your application.                   */
/*                                                                       */
/*************************************************************************/

// Include files
#include "main.h"
#include "coder_array.h"
#include "fftSqueeze.h"
#include "fftSqueeze_terminate.h"
#include "rt_nonfinite.h"

// Function Declarations
static void argInit_128x1_real_T(double result[128]);

static coder::array<double, 2U> argInit_1xUnbounded_real_T();

static double argInit_real_T();

// Function Definitions
static void argInit_128x1_real_T(double result[128])
{
  // Loop over the array to initialize each element.
  for (int idx0{0}; idx0 < 128; idx0++) {
    // Set the value of the array element.
    // Change this value to the value that the application requires.
    result[idx0] = argInit_real_T();
  }
}

static coder::array<double, 2U> argInit_1xUnbounded_real_T()
{
  coder::array<double, 2U> result;
  // Set the size of the array.
  // Change this size to the value that the application requires.
  result.set_size(1, 2);
  // Loop over the array to initialize each element.
  for (int idx0{0}; idx0 < 1; idx0++) {
    for (int idx1{0}; idx1 < result.size(1); idx1++) {
      // Set the value of the array element.
      // Change this value to the value that the application requires.
      result[idx1] = argInit_real_T();
    }
  }
  return result;
}

static double argInit_real_T()
{
  return 0.0;
}

int main(int, char **)
{
  // The initialize function is being called automatically from your entry-point
  // function. So, a call to initialize is not included here. Invoke the
  // entry-point functions.
  // You can call entry-point functions multiple times.
  main_fftSqueeze();
  // Terminate the application.
  // You do not need to do this more than one time.
  fftSqueeze_terminate();
  return 0;
}

void main_fftSqueeze()
{
  coder::array<creal_T, 2U> s;
  coder::array<double, 2U> t;
  coder::array<double, 2U> x;
  double dv[128];
  double f_data[65];
  int f_size;
  // Initialize function 'fftSqueeze' input arguments.
  // Initialize function input argument 'x'.
  x = argInit_1xUnbounded_real_T();
  // Initialize function input argument 'window'.
  // Call the entry-point 'fftSqueeze'.
  argInit_128x1_real_T(dv);
  fftSqueeze(x, argInit_real_T(), dv, s, f_data, &f_size, t);
}

// End of code generation (main.cpp)
