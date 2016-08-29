/*
 * myMetropolisHasting_initialize.c
 *
 * Code generation for function 'myMetropolisHasting_initialize'
 *
 * C source code generated on: Wed Jan 21 15:09:19 2015
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "myMetropolisHasting.h"
#include "myMetropolisHasting_initialize.h"
#include "myMetropolisHasting_data.h"

/* Function Definitions */
void myMetropolisHasting_initialize(emlrtContext *aContext)
{
  emlrtBreakCheckR2012bFlagVar = emlrtGetBreakCheckFlagAddressR2012b();
  emlrtCreateRootTLS(&emlrtRootTLSGlobal, aContext, NULL, 1);
  emlrtClearAllocCountR2012b(emlrtRootTLSGlobal, FALSE, 0U, 0);
  emlrtEnterRtStackR2012b(emlrtRootTLSGlobal);
  emlrtFirstTimeR2012b(emlrtRootTLSGlobal);
}

/* End of code generation (myMetropolisHasting_initialize.c) */
