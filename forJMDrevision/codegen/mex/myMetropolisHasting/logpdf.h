/*
 * logpdf.h
 *
 * Code generation for function 'logpdf'
 *
 * C source code generated on: Wed Jan 21 15:09:19 2015
 *
 */

#ifndef __LOGPDF_H__
#define __LOGPDF_H__
/* Include files */
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "mwmathutil.h"

#include "tmwtypes.h"
#include "mex.h"
#include "emlrt.h"
#include "blas.h"
#include "rtwtypes.h"
#include "myMetropolisHasting_types.h"

/* Function Declarations */
extern real_T logpdf(const emxArray_real_T *x, const emxArray_real_T *A, real_T C);
#ifdef __WATCOMC__
#pragma aux logpdf value [8087];
#endif
#endif
/* End of code generation (logpdf.h) */
