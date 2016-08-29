/*
 * myMetropolisHasting_mexutil.c
 *
 * Code generation for function 'myMetropolisHasting_mexutil'
 *
 * C source code generated on: Wed Jan 21 15:09:19 2015
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "myMetropolisHasting.h"
#include "myMetropolisHasting_mexutil.h"

/* Function Definitions */
const mxArray *emlrt_marshallOut(real_T u)
{
  const mxArray *y;
  const mxArray *m2;
  y = NULL;
  m2 = mxCreateDoubleScalar(u);
  emlrtAssign(&y, m2);
  return y;
}

/* End of code generation (myMetropolisHasting_mexutil.c) */
