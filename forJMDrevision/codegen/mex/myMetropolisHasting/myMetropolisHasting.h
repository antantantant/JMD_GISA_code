/*
 * myMetropolisHasting.h
 *
 * Code generation for function 'myMetropolisHasting'
 *
 * C source code generated on: Wed Jan 21 15:09:19 2015
 *
 */

#ifndef __MYMETROPOLISHASTING_H__
#define __MYMETROPOLISHASTING_H__
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
extern void eml_error(void);
extern void myMetropolisHasting(const emxArray_real_T *start, real_T nsamples, real_T burnin, const emxArray_real_T *A, real_T C, emxArray_real_T *smpl, real_T *accept);
#endif
/* End of code generation (myMetropolisHasting.h) */
