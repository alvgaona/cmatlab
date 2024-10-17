//
// Trial License - for use to evaluate programs for possible purchase as
// an end-user only.
//
// resizeim.cpp
//
// Code generation for function 'resizeim'
//

// Include files
#include "resizeim.h"
#include "imresize.h"
#include "resizeim_data.h"
#include "resizeim_initialize.h"
#include "rt_nonfinite.h"
#include "coder_array.h"

// Function Definitions
void resizeim(const coder::array<unsigned char, 3U> &im,
              coder::array<unsigned char, 3U> &out)
{
  if (!isInitialized_resizeim) {
    resizeim_initialize();
  }
  coder::imresize(im, out);
}

// End of code generation (resizeim.cpp)
