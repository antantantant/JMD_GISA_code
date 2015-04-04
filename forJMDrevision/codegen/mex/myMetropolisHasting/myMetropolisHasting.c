/*
 * myMetropolisHasting.c
 *
 * Code generation for function 'myMetropolisHasting'
 *
 * C source code generated on: Wed Jan 21 15:09:19 2015
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "myMetropolisHasting.h"
#include "myMetropolisHasting_emxutil.h"
#include "logpdf.h"
#include "myMetropolisHasting_mexutil.h"
#include "myMetropolisHasting_data.h"

/* Variable Definitions */
static emlrtRSInfo emlrtRSI = { 12, "myMetropolisHasting",
  "C:/doiUsers/Max/GGBS for DETC2014b/forJMDrevision/myMetropolisHasting.m" };

static emlrtRSInfo b_emlrtRSI = { 14, "myMetropolisHasting",
  "C:/doiUsers/Max/GGBS for DETC2014b/forJMDrevision/myMetropolisHasting.m" };

static emlrtRSInfo c_emlrtRSI = { 16, "myMetropolisHasting",
  "C:/doiUsers/Max/GGBS for DETC2014b/forJMDrevision/myMetropolisHasting.m" };

static emlrtRSInfo d_emlrtRSI = { 79, "rand",
  "C:/Program Files/MATLAB/R2013a/toolbox/eml/lib/matlab/randfun/rand.m" };

static emlrtRSInfo f_emlrtRSI = { 20, "eml_error",
  "C:/Program Files/MATLAB/R2013a/toolbox/eml/lib/matlab/eml/eml_error.m" };

static emlrtRSInfo g_emlrtRSI = { 2, "proprnd",
  "C:/doiUsers/Max/GGBS for DETC2014b/forJMDrevision/proprnd.m" };

static emlrtMCInfo emlrtMCI = { 79, 9, "rand",
  "C:/Program Files/MATLAB/R2013a/toolbox/eml/lib/matlab/randfun/rand.m" };

static emlrtRTEInfo emlrtRTEI = { 3, 26, "myMetropolisHasting",
  "C:/doiUsers/Max/GGBS for DETC2014b/forJMDrevision/myMetropolisHasting.m" };

static emlrtRTEInfo b_emlrtRTEI = { 9, 5, "myMetropolisHasting",
  "C:/doiUsers/Max/GGBS for DETC2014b/forJMDrevision/myMetropolisHasting.m" };

static emlrtRTEInfo c_emlrtRTEI = { 12, 5, "myMetropolisHasting",
  "C:/doiUsers/Max/GGBS for DETC2014b/forJMDrevision/myMetropolisHasting.m" };

static emlrtRTEInfo d_emlrtRTEI = { 14, 9, "myMetropolisHasting",
  "C:/doiUsers/Max/GGBS for DETC2014b/forJMDrevision/myMetropolisHasting.m" };

static emlrtECInfo emlrtECI = { 2, 2, 13, "proprnd",
  "C:/doiUsers/Max/GGBS for DETC2014b/forJMDrevision/proprnd.m" };

static emlrtECInfo b_emlrtECI = { -1, 25, 13, "myMetropolisHasting",
  "C:/doiUsers/Max/GGBS for DETC2014b/forJMDrevision/myMetropolisHasting.m" };

static emlrtBCInfo emlrtBCI = { -1, -1, 25, 18, "smpl", "myMetropolisHasting",
  "C:/doiUsers/Max/GGBS for DETC2014b/forJMDrevision/myMetropolisHasting.m", 0 };

static emlrtDCInfo emlrtDCI = { 25, 18, "myMetropolisHasting",
  "C:/doiUsers/Max/GGBS for DETC2014b/forJMDrevision/myMetropolisHasting.m", 1 };

static emlrtECInfo c_emlrtECI = { -1, 21, 13, "myMetropolisHasting",
  "C:/doiUsers/Max/GGBS for DETC2014b/forJMDrevision/myMetropolisHasting.m" };

static emlrtBCInfo b_emlrtBCI = { -1, -1, 18, 18, "U", "myMetropolisHasting",
  "C:/doiUsers/Max/GGBS for DETC2014b/forJMDrevision/myMetropolisHasting.m", 0 };

static emlrtDCInfo b_emlrtDCI = { 18, 18, "myMetropolisHasting",
  "C:/doiUsers/Max/GGBS for DETC2014b/forJMDrevision/myMetropolisHasting.m", 1 };

static emlrtRTEInfo g_emlrtRTEI = { 13, 5, "myMetropolisHasting",
  "C:/doiUsers/Max/GGBS for DETC2014b/forJMDrevision/myMetropolisHasting.m" };

static emlrtDCInfo c_emlrtDCI = { 7, 18, "myMetropolisHasting",
  "C:/doiUsers/Max/GGBS for DETC2014b/forJMDrevision/myMetropolisHasting.m", 1 };

static emlrtDCInfo d_emlrtDCI = { 7, 18, "myMetropolisHasting",
  "C:/doiUsers/Max/GGBS for DETC2014b/forJMDrevision/myMetropolisHasting.m", 4 };

static emlrtRTEInfo h_emlrtRTEI = { 20, 5, "eml_error",
  "C:/Program Files/MATLAB/R2013a/toolbox/eml/lib/matlab/eml/eml_error.m" };

/* Function Declarations */
static void b_emlrt_marshallIn(const mxArray *u, const emlrtMsgIdentifier
  *parentId, emxArray_real_T *y);
static const mxArray *b_rand(const mxArray *b, const mxArray *c, emlrtMCInfo
  *location);
static void emlrt_marshallIn(const mxArray *c_rand, const char_T *identifier,
  emxArray_real_T *y);
static void i_emlrt_marshallIn(const mxArray *src, const emlrtMsgIdentifier
  *msgId, emxArray_real_T *ret);

/* Function Definitions */
static void b_emlrt_marshallIn(const mxArray *u, const emlrtMsgIdentifier
  *parentId, emxArray_real_T *y)
{
  i_emlrt_marshallIn(emlrtAlias(u), parentId, y);
  emlrtDestroyArray(&u);
}

static const mxArray *b_rand(const mxArray *b, const mxArray *c, emlrtMCInfo
  *location)
{
  const mxArray *pArrays[2];
  const mxArray *m5;
  pArrays[0] = b;
  pArrays[1] = c;
  return emlrtCallMATLABR2012b(emlrtRootTLSGlobal, 1, &m5, 2, pArrays, "rand",
    TRUE, location);
}

static void emlrt_marshallIn(const mxArray *c_rand, const char_T *identifier,
  emxArray_real_T *y)
{
  emlrtMsgIdentifier thisId;
  thisId.fIdentifier = identifier;
  thisId.fParent = NULL;
  b_emlrt_marshallIn(emlrtAlias(c_rand), &thisId, y);
  emlrtDestroyArray(&c_rand);
}

static void i_emlrt_marshallIn(const mxArray *src, const emlrtMsgIdentifier
  *msgId, emxArray_real_T *ret)
{
  int32_T iv12[2];
  boolean_T bv0[2];
  int32_T i4;
  static const boolean_T bv1[2] = { FALSE, TRUE };

  int32_T iv13[2];
  for (i4 = 0; i4 < 2; i4++) {
    iv12[i4] = 1 + -2 * i4;
    bv0[i4] = bv1[i4];
  }

  emlrtCheckVsBuiltInR2012b(emlrtRootTLSGlobal, msgId, src, "double", FALSE, 2U,
    iv12, bv0, iv13);
  i4 = ret->size[0] * ret->size[1];
  ret->size[0] = iv13[0];
  ret->size[1] = iv13[1];
  emxEnsureCapacity((emxArray__common *)ret, i4, (int32_T)sizeof(real_T),
                    (emlrtRTEInfo *)NULL);
  emlrtImportArrayR2011b(src, ret->data, 8, FALSE);
  emlrtDestroyArray(&src);
}

void eml_error(void)
{
  static char_T cv0[3][1] = { { 'l' }, { 'o' }, { 'g' } };

  emlrtPushRtStackR2012b(&f_emlrtRSI, emlrtRootTLSGlobal);
  emlrtErrorWithMessageIdR2012b(emlrtRootTLSGlobal, &h_emlrtRTEI,
    "Coder:toolbox:ElFunDomainError", 3, 4, 3, cv0);
  emlrtPopRtStackR2012b(&f_emlrtRSI, emlrtRootTLSGlobal);
}

void myMetropolisHasting(const emxArray_real_T *start, real_T nsamples, real_T
  burnin, const emxArray_real_T *A, real_T C, emxArray_real_T *smpl, real_T
  *accept)
{
  real_T b_nsamples[2];
  real_T dv0[2];
  int32_T i0;
  real_T d0;
  int32_T loop_ub;
  emxArray_real_T *x0;
  emxArray_real_T *U;
  emxArray_real_T *y;
  int32_T i;
  emxArray_int32_T *r0;
  emxArray_real_T *b_y;
  emxArray_real_T *b_x0;
  emxArray_real_T *c_y;
  real_T b_i;
  int32_T i1;
  int32_T acc;
  int32_T c_x0[2];
  int32_T d_y[2];
  real_T rho;
  int32_T iv0[2];
  int32_T iv1[2];
  emlrtHeapReferenceStackEnterFcnR2012b(emlrtRootTLSGlobal);

  /* my metropolis hasting algorithm */
  /*  Put the replicates dimension second. */
  b_nsamples[0] = nsamples;
  b_nsamples[1] = start->size[1];
  for (i0 = 0; i0 < 2; i0++) {
    d0 = emlrtNonNegativeCheckFastR2012b(b_nsamples[i0], &d_emlrtDCI,
      emlrtRootTLSGlobal);
    dv0[i0] = emlrtIntegerCheckFastR2012b(d0, &c_emlrtDCI, emlrtRootTLSGlobal);
  }

  i0 = smpl->size[0] * smpl->size[1];
  smpl->size[0] = (int32_T)dv0[0];
  emxEnsureCapacity((emxArray__common *)smpl, i0, (int32_T)sizeof(real_T),
                    &emlrtRTEI);
  i0 = smpl->size[0] * smpl->size[1];
  smpl->size[1] = (int32_T)dv0[1];
  emxEnsureCapacity((emxArray__common *)smpl, i0, (int32_T)sizeof(real_T),
                    &emlrtRTEI);
  loop_ub = (int32_T)dv0[0] * (int32_T)dv0[1];
  for (i0 = 0; i0 < loop_ub; i0++) {
    smpl->data[i0] = 0.0;
  }

  emxInit_real_T(&x0, 2, &b_emlrtRTEI, TRUE);
  i0 = x0->size[0] * x0->size[1];
  x0->size[0] = 1;
  x0->size[1] = start->size[1];
  emxEnsureCapacity((emxArray__common *)x0, i0, (int32_T)sizeof(real_T),
                    &emlrtRTEI);
  loop_ub = start->size[0] * start->size[1];
  for (i0 = 0; i0 < loop_ub; i0++) {
    x0->data[i0] = start->data[i0];
  }

  emxInit_real_T(&U, 2, &c_emlrtRTEI, TRUE);
  emxInit_real_T(&y, 2, &d_emlrtRTEI, TRUE);

  /* x0  is the place holder for the current value */
  *accept = 0.0;

  /*  Metropolis-Hasting Algorithm. */
  emlrtPushRtStackR2012b(&emlrtRSI, emlrtRootTLSGlobal);
  emlrtPushRtStackR2012b(&d_emlrtRSI, emlrtRootTLSGlobal);
  emlrt_marshallIn(b_rand(emlrt_marshallOut(1.0), emlrt_marshallOut(nsamples +
    burnin), &emlrtMCI), "rand", y);
  emlrtPopRtStackR2012b(&d_emlrtRSI, emlrtRootTLSGlobal);
  i0 = U->size[0] * U->size[1];
  U->size[0] = 1;
  U->size[1] = y->size[1];
  emxEnsureCapacity((emxArray__common *)U, i0, (int32_T)sizeof(real_T),
                    &emlrtRTEI);
  loop_ub = y->size[0] * y->size[1];
  for (i0 = 0; i0 < loop_ub; i0++) {
    U->data[i0] = y->data[i0];
  }

  for (loop_ub = 0; loop_ub < y->size[1]; loop_ub++) {
    if (y->data[(int32_T)(1.0 + (real_T)loop_ub) - 1] < 0.0) {
      emlrtPushRtStackR2012b(&e_emlrtRSI, emlrtRootTLSGlobal);
      eml_error();
      emlrtPopRtStackR2012b(&e_emlrtRSI, emlrtRootTLSGlobal);
    }
  }

  for (loop_ub = 0; loop_ub < y->size[1]; loop_ub++) {
    U->data[(int32_T)(1.0 + (real_T)loop_ub) - 1] = muDoubleScalarLog(U->data
      [(int32_T)(1.0 + (real_T)loop_ub) - 1]);
  }

  emlrtPopRtStackR2012b(&emlrtRSI, emlrtRootTLSGlobal);
  i0 = (int32_T)(nsamples + (1.0 - (1.0 - burnin)));
  emlrtForLoopVectorCheckR2012b(1.0 - burnin, 1.0, nsamples, mxDOUBLE_CLASS, i0,
    &g_emlrtRTEI, emlrtRootTLSGlobal);
  i = 0;
  emxInit_int32_T(&r0, 1, &emlrtRTEI, TRUE);
  b_emxInit_real_T(&b_y, 1, &emlrtRTEI, TRUE);
  b_emxInit_real_T(&b_x0, 1, &emlrtRTEI, TRUE);
  emxInit_real_T(&c_y, 2, &emlrtRTEI, TRUE);
  while (i <= i0 - 1) {
    b_i = (1.0 - burnin) + (real_T)i;
    emlrtPushRtStackR2012b(&b_emlrtRSI, emlrtRootTLSGlobal);
    emlrtPushRtStackR2012b(&g_emlrtRSI, emlrtRootTLSGlobal);
    emlrtPushRtStackR2012b(&d_emlrtRSI, emlrtRootTLSGlobal);
    emlrt_marshallIn(b_rand(emlrt_marshallOut(1.0), emlrt_marshallOut(x0->size[1]),
      &emlrtMCI), "rand", y);
    emlrtPopRtStackR2012b(&d_emlrtRSI, emlrtRootTLSGlobal);
    i1 = y->size[0] * y->size[1];
    y->size[0] = 1;
    emxEnsureCapacity((emxArray__common *)y, i1, (int32_T)sizeof(real_T),
                      &emlrtRTEI);
    loop_ub = y->size[0];
    acc = y->size[1];
    loop_ub *= acc;
    for (i1 = 0; i1 < loop_ub; i1++) {
      y->data[i1] = (y->data[i1] - 0.5) * 0.1;
    }

    emlrtPopRtStackR2012b(&g_emlrtRSI, emlrtRootTLSGlobal);
    for (i1 = 0; i1 < 2; i1++) {
      c_x0[i1] = x0->size[i1];
    }

    for (i1 = 0; i1 < 2; i1++) {
      d_y[i1] = y->size[i1];
    }

    emlrtSizeEqCheck2DFastR2012b(c_x0, d_y, &emlrtECI, emlrtRootTLSGlobal);
    i1 = y->size[0] * y->size[1];
    y->size[0] = 1;
    y->size[1] = x0->size[1];
    emxEnsureCapacity((emxArray__common *)y, i1, (int32_T)sizeof(real_T),
                      &emlrtRTEI);
    loop_ub = x0->size[0] * x0->size[1];
    for (i1 = 0; i1 < loop_ub; i1++) {
      y->data[i1] += x0->data[i1];
    }

    emlrtPopRtStackR2012b(&b_emlrtRSI, emlrtRootTLSGlobal);

    /*  sample from proposal dist'n */
    /* save the evaluation time for symmetric proposal dist'n */
    emlrtPushRtStackR2012b(&c_emlrtRSI, emlrtRootTLSGlobal);
    i1 = b_y->size[0];
    b_y->size[0] = y->size[1];
    emxEnsureCapacity((emxArray__common *)b_y, i1, (int32_T)sizeof(real_T),
                      &emlrtRTEI);
    loop_ub = y->size[1];
    for (i1 = 0; i1 < loop_ub; i1++) {
      b_y->data[i1] = y->data[i1];
    }

    i1 = b_x0->size[0];
    b_x0->size[0] = x0->size[1];
    emxEnsureCapacity((emxArray__common *)b_x0, i1, (int32_T)sizeof(real_T),
                      &emlrtRTEI);
    loop_ub = x0->size[1];
    for (i1 = 0; i1 < loop_ub; i1++) {
      b_x0->data[i1] = x0->data[i1];
    }

    rho = logpdf(b_y, A, C) - logpdf(b_x0, A, C);
    emlrtPopRtStackR2012b(&c_emlrtRSI, emlrtRootTLSGlobal);

    /*  Accept or reject the proposal. */
    i1 = U->size[1];
    d0 = b_i + burnin;
    loop_ub = (int32_T)emlrtIntegerCheckFastR2012b(d0, &b_emlrtDCI,
      emlrtRootTLSGlobal);
    emlrtDynamicBoundsCheckFastR2012b(loop_ub, 1, i1, &b_emlrtBCI,
      emlrtRootTLSGlobal);
    acc = (U->data[U->size[0] * ((int32_T)(b_i + burnin) - 1)] <=
           muDoubleScalarMin(rho, 0.0));
    if (acc > 0) {
      loop_ub = x0->size[1];
      i1 = r0->size[0];
      r0->size[0] = loop_ub;
      emxEnsureCapacity((emxArray__common *)r0, i1, (int32_T)sizeof(int32_T),
                        &emlrtRTEI);
      for (i1 = 0; i1 < loop_ub; i1++) {
        r0->data[i1] = i1;
      }

      iv0[0] = 1;
      iv0[1] = r0->size[0];
      loop_ub = y->size[1];
      i1 = c_y->size[0] * c_y->size[1];
      c_y->size[0] = 1;
      c_y->size[1] = loop_ub;
      emxEnsureCapacity((emxArray__common *)c_y, i1, (int32_T)sizeof(real_T),
                        &emlrtRTEI);
      for (i1 = 0; i1 < loop_ub; i1++) {
        c_y->data[c_y->size[0] * i1] = y->data[y->size[0] * i1];
      }

      for (i1 = 0; i1 < 2; i1++) {
        d_y[i1] = c_y->size[i1];
      }

      emlrtSubAssignSizeCheckR2012b(iv0, 2, d_y, 2, &c_emlrtECI,
        emlrtRootTLSGlobal);
      loop_ub = y->size[1] - 1;
      for (i1 = 0; i1 <= loop_ub; i1++) {
        x0->data[x0->size[0] * r0->data[i1]] = y->data[y->size[0] * i1];
      }

      /*  preserves x's shape. */
    }

    *accept += (real_T)acc;
    if (b_i > 0.0) {
      /*  burnin */
      i1 = smpl->size[0];
      loop_ub = (int32_T)emlrtIntegerCheckFastR2012b(b_i, &emlrtDCI,
        emlrtRootTLSGlobal);
      emlrtDynamicBoundsCheckFastR2012b(loop_ub, 1, i1, &emlrtBCI,
        emlrtRootTLSGlobal);
      loop_ub = smpl->size[1];
      i1 = r0->size[0];
      r0->size[0] = loop_ub;
      emxEnsureCapacity((emxArray__common *)r0, i1, (int32_T)sizeof(int32_T),
                        &emlrtRTEI);
      for (i1 = 0; i1 < loop_ub; i1++) {
        r0->data[i1] = i1;
      }

      iv1[0] = 1;
      iv1[1] = r0->size[0];
      emlrtSubAssignSizeCheckR2012b(iv1, 2, *(int32_T (*)[2])x0->size, 2,
        &b_emlrtECI, emlrtRootTLSGlobal);
      loop_ub = x0->size[1];
      for (i1 = 0; i1 < loop_ub; i1++) {
        smpl->data[((int32_T)b_i + smpl->size[0] * r0->data[i1]) - 1] = x0->
          data[x0->size[0] * i1];
      }
    }

    i++;
    emlrtBreakCheckFastR2012b(emlrtBreakCheckR2012bFlagVar, emlrtRootTLSGlobal);
  }

  emxFree_real_T(&c_y);
  emxFree_real_T(&b_x0);
  emxFree_real_T(&b_y);
  emxFree_int32_T(&r0);
  emxFree_real_T(&y);
  emxFree_real_T(&U);
  emxFree_real_T(&x0);

  /*  Accept rate can be used to optimize the choice of scale parameters in */
  /*  random walk MH sampler. See for example Roberts, Gelman and Gilks (1997). */
  *accept /= nsamples + burnin;
  emlrtHeapReferenceStackLeaveFcnR2012b(emlrtRootTLSGlobal);
}

/* End of code generation (myMetropolisHasting.c) */
