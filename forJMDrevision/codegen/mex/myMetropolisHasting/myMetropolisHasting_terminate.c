/*
 * myMetropolisHasting_terminate.c
 *
 * Code generation for function 'myMetropolisHasting_terminate'
 *
 * C source code generated on: Wed Jan 21 15:09:19 2015
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "myMetropolisHasting.h"
#include "myMetropolisHasting_terminate.h"

/* Function Definitions */
void myMetropolisHasting_atexit(void)
{
  emlrtCreateRootTLS(&emlrtRootTLSGlobal, &emlrtContextGlobal, NULL, 1);
  emlrtEnterRtStackR2012b(emlrtRootTLSGlobal);
  emlrtLeaveRtStackR2012b(emlrtRootTLSGlobal);
  emlrtDestroyRootTLS(&emlrtRootTLSGlobal);
}

void myMetropolisHasting_terminate(void)
{
  emlrtLeaveRtStackR2012b(emlrtRootTLSGlobal);
  emlrtDestroyRootTLS(&emlrtRootTLSGlobal);
}

/* End of code generation (myMetropolisHasting_terminate.c) */
