/*
 * logpdf.c
 *
 * Code generation for function 'logpdf'
 *
 * C source code generated on: Wed Jan 21 15:09:19 2015
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "myMetropolisHasting.h"
#include "logpdf.h"
#include "myMetropolisHasting_emxutil.h"
#include "myMetropolisHasting_data.h"

/* Variable Definitions */
static emlrtRSInfo h_emlrtRSI = { 2, "logpdf",
  "C:/doiUsers/Max/GGBS for DETC2014b/forJMDrevision/logpdf.m" };

static emlrtRSInfo i_emlrtRSI = { 55, "mtimes",
  "C:/Program Files/MATLAB/R2013a/toolbox/eml/lib/matlab/ops/mtimes.m" };

static emlrtRSInfo j_emlrtRSI = { 21, "mtimes",
  "C:/Program Files/MATLAB/R2013a/toolbox/eml/lib/matlab/ops/mtimes.m" };

static emlrtRSInfo k_emlrtRSI = { 89, "mtimes",
  "C:/Program Files/MATLAB/R2013a/toolbox/eml/lib/matlab/ops/mtimes.m" };

static emlrtRSInfo l_emlrtRSI = { 84, "mtimes",
  "C:/Program Files/MATLAB/R2013a/toolbox/eml/lib/matlab/ops/mtimes.m" };

static emlrtRSInfo m_emlrtRSI = { 54, "eml_xgemm",
  "C:/Program Files/MATLAB/R2013a/toolbox/eml/lib/matlab/eml/blas/eml_xgemm.m" };

static emlrtRSInfo o_emlrtRSI = { 32, "eml_blas_xgemm",
  "C:/Program Files/MATLAB/R2013a/toolbox/eml/lib/matlab/eml/blas/external/eml_blas_xgemm.m"
};

static emlrtRSInfo s_emlrtRSI = { 51, "eml_int_forloop_overflow_check",
  "C:/Program Files/MATLAB/R2013a/toolbox/eml/lib/matlab/eml/eml_int_forloop_overflow_check.m"
};

static emlrtRSInfo t_emlrtRSI = { 12, "eml_int_forloop_overflow_check",
  "C:/Program Files/MATLAB/R2013a/toolbox/eml/lib/matlab/eml/eml_int_forloop_overflow_check.m"
};

static emlrtRSInfo u_emlrtRSI = { 110, "eml_blas_xgemm",
  "C:/Program Files/MATLAB/R2013a/toolbox/eml/lib/matlab/eml/blas/external/eml_blas_xgemm.m"
};

static emlrtRSInfo v_emlrtRSI = { 111, "eml_blas_xgemm",
  "C:/Program Files/MATLAB/R2013a/toolbox/eml/lib/matlab/eml/blas/external/eml_blas_xgemm.m"
};

static emlrtRSInfo w_emlrtRSI = { 112, "eml_blas_xgemm",
  "C:/Program Files/MATLAB/R2013a/toolbox/eml/lib/matlab/eml/blas/external/eml_blas_xgemm.m"
};

static emlrtRSInfo x_emlrtRSI = { 113, "eml_blas_xgemm",
  "C:/Program Files/MATLAB/R2013a/toolbox/eml/lib/matlab/eml/blas/external/eml_blas_xgemm.m"
};

static emlrtRSInfo y_emlrtRSI = { 114, "eml_blas_xgemm",
  "C:/Program Files/MATLAB/R2013a/toolbox/eml/lib/matlab/eml/blas/external/eml_blas_xgemm.m"
};

static emlrtRSInfo ab_emlrtRSI = { 115, "eml_blas_xgemm",
  "C:/Program Files/MATLAB/R2013a/toolbox/eml/lib/matlab/eml/blas/external/eml_blas_xgemm.m"
};

static emlrtRSInfo bb_emlrtRSI = { 119, "eml_blas_xgemm",
  "C:/Program Files/MATLAB/R2013a/toolbox/eml/lib/matlab/eml/blas/external/eml_blas_xgemm.m"
};

static emlrtRSInfo cb_emlrtRSI = { 122, "eml_blas_xgemm",
  "C:/Program Files/MATLAB/R2013a/toolbox/eml/lib/matlab/eml/blas/external/eml_blas_xgemm.m"
};

static emlrtRSInfo db_emlrtRSI = { 125, "eml_blas_xgemm",
  "C:/Program Files/MATLAB/R2013a/toolbox/eml/lib/matlab/eml/blas/external/eml_blas_xgemm.m"
};

static emlrtRSInfo eb_emlrtRSI = { 128, "eml_blas_xgemm",
  "C:/Program Files/MATLAB/R2013a/toolbox/eml/lib/matlab/eml/blas/external/eml_blas_xgemm.m"
};

static emlrtRSInfo fb_emlrtRSI = { 131, "eml_blas_xgemm",
  "C:/Program Files/MATLAB/R2013a/toolbox/eml/lib/matlab/eml/blas/external/eml_blas_xgemm.m"
};

static emlrtRSInfo gb_emlrtRSI = { 134, "eml_blas_xgemm",
  "C:/Program Files/MATLAB/R2013a/toolbox/eml/lib/matlab/eml/blas/external/eml_blas_xgemm.m"
};

static emlrtRSInfo hb_emlrtRSI = { 14, "eml_c_cast",
  "C:/Program Files/MATLAB/R2013a/toolbox/eml/lib/matlab/eml/blas/external/eml_c_cast.m"
};

static emlrtRSInfo ib_emlrtRSI = { 17, "sum",
  "C:/Program Files/MATLAB/R2013a/toolbox/eml/lib/matlab/datafun/sum.m" };

static emlrtRSInfo jb_emlrtRSI = { 20, "sum",
  "C:/Program Files/MATLAB/R2013a/toolbox/eml/lib/matlab/datafun/sum.m" };

static emlrtRSInfo kb_emlrtRSI = { 61, "sum",
  "C:/Program Files/MATLAB/R2013a/toolbox/eml/lib/matlab/datafun/sum.m" };

static emlrtRSInfo lb_emlrtRSI = { 49, "mtimes",
  "C:/Program Files/MATLAB/R2013a/toolbox/eml/lib/matlab/ops/mtimes.m" };

static emlrtRSInfo mb_emlrtRSI = { 31, "eml_xdotu",
  "C:/Program Files/MATLAB/R2013a/toolbox/eml/lib/matlab/eml/blas/eml_xdotu.m" };

static emlrtRSInfo nb_emlrtRSI = { 28, "eml_xdot",
  "C:/Program Files/MATLAB/R2013a/toolbox/eml/lib/matlab/eml/blas/eml_xdot.m" };

static emlrtRSInfo pb_emlrtRSI = { 28, "eml_blas_xdot",
  "C:/Program Files/MATLAB/R2013a/toolbox/eml/lib/matlab/eml/blas/external/eml_blas_xdot.m"
};

static emlrtRSInfo sb_emlrtRSI = { 64, "eml_blas_xdot",
  "C:/Program Files/MATLAB/R2013a/toolbox/eml/lib/matlab/eml/blas/external/eml_blas_xdot.m"
};

static emlrtRSInfo tb_emlrtRSI = { 65, "eml_blas_xdot",
  "C:/Program Files/MATLAB/R2013a/toolbox/eml/lib/matlab/eml/blas/external/eml_blas_xdot.m"
};

static emlrtRSInfo ub_emlrtRSI = { 66, "eml_blas_xdot",
  "C:/Program Files/MATLAB/R2013a/toolbox/eml/lib/matlab/eml/blas/external/eml_blas_xdot.m"
};

static emlrtRSInfo vb_emlrtRSI = { 70, "eml_blas_xdot",
  "C:/Program Files/MATLAB/R2013a/toolbox/eml/lib/matlab/eml/blas/external/eml_blas_xdot.m"
};

static emlrtRSInfo wb_emlrtRSI = { 73, "eml_blas_xdot",
  "C:/Program Files/MATLAB/R2013a/toolbox/eml/lib/matlab/eml/blas/external/eml_blas_xdot.m"
};

static emlrtMCInfo b_emlrtMCI = { 85, 13, "mtimes",
  "C:/Program Files/MATLAB/R2013a/toolbox/eml/lib/matlab/ops/mtimes.m" };

static emlrtMCInfo c_emlrtMCI = { 84, 23, "mtimes",
  "C:/Program Files/MATLAB/R2013a/toolbox/eml/lib/matlab/ops/mtimes.m" };

static emlrtMCInfo d_emlrtMCI = { 90, 13, "mtimes",
  "C:/Program Files/MATLAB/R2013a/toolbox/eml/lib/matlab/ops/mtimes.m" };

static emlrtMCInfo e_emlrtMCI = { 89, 23, "mtimes",
  "C:/Program Files/MATLAB/R2013a/toolbox/eml/lib/matlab/ops/mtimes.m" };

static emlrtMCInfo f_emlrtMCI = { 52, 9, "eml_int_forloop_overflow_check",
  "C:/Program Files/MATLAB/R2013a/toolbox/eml/lib/matlab/eml/eml_int_forloop_overflow_check.m"
};

static emlrtMCInfo g_emlrtMCI = { 51, 15, "eml_int_forloop_overflow_check",
  "C:/Program Files/MATLAB/R2013a/toolbox/eml/lib/matlab/eml/eml_int_forloop_overflow_check.m"
};

static emlrtMCInfo h_emlrtMCI = { 18, 9, "sum",
  "C:/Program Files/MATLAB/R2013a/toolbox/eml/lib/matlab/datafun/sum.m" };

static emlrtMCInfo i_emlrtMCI = { 17, 19, "sum",
  "C:/Program Files/MATLAB/R2013a/toolbox/eml/lib/matlab/datafun/sum.m" };

static emlrtMCInfo j_emlrtMCI = { 23, 9, "sum",
  "C:/Program Files/MATLAB/R2013a/toolbox/eml/lib/matlab/datafun/sum.m" };

static emlrtMCInfo k_emlrtMCI = { 20, 19, "sum",
  "C:/Program Files/MATLAB/R2013a/toolbox/eml/lib/matlab/datafun/sum.m" };

static emlrtRTEInfo e_emlrtRTEI = { 1, 14, "logpdf",
  "C:/doiUsers/Max/GGBS for DETC2014b/forJMDrevision/logpdf.m" };

/* Function Declarations */
static const mxArray *b_message(const mxArray *b, const mxArray *c, emlrtMCInfo *
  location);
static void check_forloop_overflow_error(void);
static void error(const mxArray *b, emlrtMCInfo *location);
static const mxArray *message(const mxArray *b, emlrtMCInfo *location);

/* Function Definitions */
static const mxArray *b_message(const mxArray *b, const mxArray *c, emlrtMCInfo *
  location)
{
  const mxArray *pArrays[2];
  const mxArray *m7;
  pArrays[0] = b;
  pArrays[1] = c;
  return emlrtCallMATLABR2012b(emlrtRootTLSGlobal, 1, &m7, 2, pArrays, "message",
    TRUE, location);
}

static void check_forloop_overflow_error(void)
{
  const mxArray *y;
  static const int32_T iv8[2] = { 1, 34 };

  const mxArray *m1;
  char_T cv9[34];
  int32_T i;
  static const char_T cv10[34] = { 'C', 'o', 'd', 'e', 'r', ':', 't', 'o', 'o',
    'l', 'b', 'o', 'x', ':', 'i', 'n', 't', '_', 'f', 'o', 'r', 'l', 'o', 'o',
    'p', '_', 'o', 'v', 'e', 'r', 'f', 'l', 'o', 'w' };

  const mxArray *b_y;
  static const int32_T iv9[2] = { 1, 23 };

  char_T cv11[23];
  static const char_T cv12[23] = { 'c', 'o', 'd', 'e', 'r', '.', 'i', 'n', 't',
    'e', 'r', 'n', 'a', 'l', '.', 'i', 'n', 'd', 'e', 'x', 'I', 'n', 't' };

  emlrtPushRtStackR2012b(&s_emlrtRSI, emlrtRootTLSGlobal);
  y = NULL;
  m1 = mxCreateCharArray(2, iv8);
  for (i = 0; i < 34; i++) {
    cv9[i] = cv10[i];
  }

  emlrtInitCharArrayR2013a(emlrtRootTLSGlobal, 34, m1, cv9);
  emlrtAssign(&y, m1);
  b_y = NULL;
  m1 = mxCreateCharArray(2, iv9);
  for (i = 0; i < 23; i++) {
    cv11[i] = cv12[i];
  }

  emlrtInitCharArrayR2013a(emlrtRootTLSGlobal, 23, m1, cv11);
  emlrtAssign(&b_y, m1);
  error(b_message(y, b_y, &f_emlrtMCI), &g_emlrtMCI);
  emlrtPopRtStackR2012b(&s_emlrtRSI, emlrtRootTLSGlobal);
}

static void error(const mxArray *b, emlrtMCInfo *location)
{
  const mxArray *pArray;
  pArray = b;
  emlrtCallMATLABR2012b(emlrtRootTLSGlobal, 0, NULL, 1, &pArray, "error", TRUE,
                        location);
}

static const mxArray *message(const mxArray *b, emlrtMCInfo *location)
{
  const mxArray *pArray;
  const mxArray *m6;
  pArray = b;
  return emlrtCallMATLABR2012b(emlrtRootTLSGlobal, 1, &m6, 1, &pArray, "message",
    TRUE, location);
}

real_T logpdf(const emxArray_real_T *x, const emxArray_real_T *A, real_T C)
{
  real_T f;
  emxArray_real_T *a;
  int32_T i2;
  int32_T i;
  const mxArray *y;
  static const int32_T iv2[2] = { 1, 45 };

  const mxArray *m0;
  char_T cv1[45];
  static const char_T cv2[45] = { 'C', 'o', 'd', 'e', 'r', ':', 't', 'o', 'o',
    'l', 'b', 'o', 'x', ':', 'm', 't', 'i', 'm', 'e', 's', '_', 'n', 'o', 'D',
    'y', 'n', 'a', 'm', 'i', 'c', 'S', 'c', 'a', 'l', 'a', 'r', 'E', 'x', 'p',
    'a', 'n', 's', 'i', 'o', 'n' };

  const mxArray *b_y;
  static const int32_T iv3[2] = { 1, 21 };

  char_T cv3[21];
  static const char_T cv4[21] = { 'C', 'o', 'd', 'e', 'r', ':', 'M', 'A', 'T',
    'L', 'A', 'B', ':', 'i', 'n', 'n', 'e', 'r', 'd', 'i', 'm' };

  emxArray_real_T *c_y;
  int32_T loop_ub;
  int32_T i3;
  uint32_T unnamed_idx_0;
  real_T alpha1;
  real_T beta1;
  char_T TRANSB;
  char_T TRANSA;
  ptrdiff_t m_t;
  ptrdiff_t n_t;
  ptrdiff_t k_t;
  ptrdiff_t lda_t;
  ptrdiff_t ldb_t;
  ptrdiff_t ldc_t;
  double * alpha1_t;
  double * Aia0_t;
  double * Bib0_t;
  double * beta1_t;
  double * Cic0_t;
  emxArray_real_T *b_x;
  boolean_T overflow;
  boolean_T p;
  int32_T exitg1;
  const mxArray *d_y;
  static const int32_T iv4[2] = { 1, 30 };

  char_T cv5[30];
  static const char_T cv6[30] = { 'C', 'o', 'd', 'e', 'r', ':', 't', 'o', 'o',
    'l', 'b', 'o', 'x', ':', 's', 'u', 'm', '_', 's', 'p', 'e', 'c', 'i', 'a',
    'l', 'E', 'm', 'p', 't', 'y' };

  const mxArray *e_y;
  static const int32_T iv5[2] = { 1, 36 };

  char_T cv7[36];
  static const char_T cv8[36] = { 'C', 'o', 'd', 'e', 'r', ':', 't', 'o', 'o',
    'l', 'b', 'o', 'x', ':', 'a', 'u', 't', 'o', 'D', 'i', 'm', 'I', 'n', 'c',
    'o', 'm', 'p', 'a', 't', 'i', 'b', 'i', 'l', 'i', 't', 'y' };

  emxArray_real_T *b_a;
  const mxArray *f_y;
  static const int32_T iv6[2] = { 1, 45 };

  const mxArray *g_y;
  static const int32_T iv7[2] = { 1, 21 };

  emlrtHeapReferenceStackEnterFcnR2012b(emlrtRootTLSGlobal);
  emxInit_real_T(&a, 2, &e_emlrtRTEI, TRUE);
  emlrtPushRtStackR2012b(&h_emlrtRSI, emlrtRootTLSGlobal);
  i2 = a->size[0] * a->size[1];
  a->size[0] = A->size[0];
  a->size[1] = A->size[1];
  emxEnsureCapacity((emxArray__common *)a, i2, (int32_T)sizeof(real_T),
                    &e_emlrtRTEI);
  i = A->size[0] * A->size[1];
  for (i2 = 0; i2 < i; i2++) {
    a->data[i2] = -A->data[i2];
  }

  emlrtPushRtStackR2012b(&j_emlrtRSI, emlrtRootTLSGlobal);
  if (!(a->size[1] == x->size[0])) {
    if (((a->size[0] == 1) && (a->size[1] == 1)) || (x->size[0] == 1)) {
      emlrtPushRtStackR2012b(&l_emlrtRSI, emlrtRootTLSGlobal);
      y = NULL;
      m0 = mxCreateCharArray(2, iv2);
      for (i = 0; i < 45; i++) {
        cv1[i] = cv2[i];
      }

      emlrtInitCharArrayR2013a(emlrtRootTLSGlobal, 45, m0, cv1);
      emlrtAssign(&y, m0);
      error(message(y, &b_emlrtMCI), &c_emlrtMCI);
      emlrtPopRtStackR2012b(&l_emlrtRSI, emlrtRootTLSGlobal);
    } else {
      emlrtPushRtStackR2012b(&k_emlrtRSI, emlrtRootTLSGlobal);
      b_y = NULL;
      m0 = mxCreateCharArray(2, iv3);
      for (i = 0; i < 21; i++) {
        cv3[i] = cv4[i];
      }

      emlrtInitCharArrayR2013a(emlrtRootTLSGlobal, 21, m0, cv3);
      emlrtAssign(&b_y, m0);
      error(message(b_y, &d_emlrtMCI), &e_emlrtMCI);
      emlrtPopRtStackR2012b(&k_emlrtRSI, emlrtRootTLSGlobal);
    }
  }

  emlrtPopRtStackR2012b(&j_emlrtRSI, emlrtRootTLSGlobal);
  b_emxInit_real_T(&c_y, 1, &e_emlrtRTEI, TRUE);
  if ((a->size[1] == 1) || (x->size[0] == 1)) {
    i2 = c_y->size[0];
    c_y->size[0] = a->size[0];
    emxEnsureCapacity((emxArray__common *)c_y, i2, (int32_T)sizeof(real_T),
                      &e_emlrtRTEI);
    i = a->size[0];
    for (i2 = 0; i2 < i; i2++) {
      c_y->data[i2] = 0.0;
      loop_ub = a->size[1];
      for (i3 = 0; i3 < loop_ub; i3++) {
        c_y->data[i2] += a->data[i2 + a->size[0] * i3] * x->data[i3];
      }
    }
  } else {
    unnamed_idx_0 = (uint32_T)a->size[0];
    emlrtPushRtStackR2012b(&i_emlrtRSI, emlrtRootTLSGlobal);
    emlrtPushRtStackR2012b(&m_emlrtRSI, emlrtRootTLSGlobal);
    i2 = c_y->size[0];
    c_y->size[0] = (int32_T)unnamed_idx_0;
    emxEnsureCapacity((emxArray__common *)c_y, i2, (int32_T)sizeof(real_T),
                      &e_emlrtRTEI);
    i = (int32_T)unnamed_idx_0;
    for (i2 = 0; i2 < i; i2++) {
      c_y->data[i2] = 0.0;
    }

    if ((a->size[0] < 1) || (a->size[1] < 1)) {
    } else {
      emlrtPushRtStackR2012b(&o_emlrtRSI, emlrtRootTLSGlobal);
      alpha1 = 1.0;
      beta1 = 0.0;
      TRANSB = 'N';
      TRANSA = 'N';
      emlrtPushRtStackR2012b(&u_emlrtRSI, emlrtRootTLSGlobal);
      emlrtPushRtStackR2012b(&hb_emlrtRSI, emlrtRootTLSGlobal);
      m_t = (ptrdiff_t)(a->size[0]);
      emlrtPopRtStackR2012b(&hb_emlrtRSI, emlrtRootTLSGlobal);
      emlrtPopRtStackR2012b(&u_emlrtRSI, emlrtRootTLSGlobal);
      emlrtPushRtStackR2012b(&v_emlrtRSI, emlrtRootTLSGlobal);
      emlrtPushRtStackR2012b(&hb_emlrtRSI, emlrtRootTLSGlobal);
      n_t = (ptrdiff_t)(1);
      emlrtPopRtStackR2012b(&hb_emlrtRSI, emlrtRootTLSGlobal);
      emlrtPopRtStackR2012b(&v_emlrtRSI, emlrtRootTLSGlobal);
      emlrtPushRtStackR2012b(&w_emlrtRSI, emlrtRootTLSGlobal);
      emlrtPushRtStackR2012b(&hb_emlrtRSI, emlrtRootTLSGlobal);
      k_t = (ptrdiff_t)(a->size[1]);
      emlrtPopRtStackR2012b(&hb_emlrtRSI, emlrtRootTLSGlobal);
      emlrtPopRtStackR2012b(&w_emlrtRSI, emlrtRootTLSGlobal);
      emlrtPushRtStackR2012b(&x_emlrtRSI, emlrtRootTLSGlobal);
      emlrtPushRtStackR2012b(&hb_emlrtRSI, emlrtRootTLSGlobal);
      lda_t = (ptrdiff_t)(a->size[0]);
      emlrtPopRtStackR2012b(&hb_emlrtRSI, emlrtRootTLSGlobal);
      emlrtPopRtStackR2012b(&x_emlrtRSI, emlrtRootTLSGlobal);
      emlrtPushRtStackR2012b(&y_emlrtRSI, emlrtRootTLSGlobal);
      emlrtPushRtStackR2012b(&hb_emlrtRSI, emlrtRootTLSGlobal);
      ldb_t = (ptrdiff_t)(a->size[1]);
      emlrtPopRtStackR2012b(&hb_emlrtRSI, emlrtRootTLSGlobal);
      emlrtPopRtStackR2012b(&y_emlrtRSI, emlrtRootTLSGlobal);
      emlrtPushRtStackR2012b(&ab_emlrtRSI, emlrtRootTLSGlobal);
      emlrtPushRtStackR2012b(&hb_emlrtRSI, emlrtRootTLSGlobal);
      ldc_t = (ptrdiff_t)(a->size[0]);
      emlrtPopRtStackR2012b(&hb_emlrtRSI, emlrtRootTLSGlobal);
      emlrtPopRtStackR2012b(&ab_emlrtRSI, emlrtRootTLSGlobal);
      emlrtPushRtStackR2012b(&bb_emlrtRSI, emlrtRootTLSGlobal);
      alpha1_t = (double *)(&alpha1);
      emlrtPopRtStackR2012b(&bb_emlrtRSI, emlrtRootTLSGlobal);
      emlrtPushRtStackR2012b(&cb_emlrtRSI, emlrtRootTLSGlobal);
      Aia0_t = (double *)(&a->data[0]);
      emlrtPopRtStackR2012b(&cb_emlrtRSI, emlrtRootTLSGlobal);
      emlrtPushRtStackR2012b(&db_emlrtRSI, emlrtRootTLSGlobal);
      Bib0_t = (double *)(&x->data[0]);
      emlrtPopRtStackR2012b(&db_emlrtRSI, emlrtRootTLSGlobal);
      emlrtPushRtStackR2012b(&eb_emlrtRSI, emlrtRootTLSGlobal);
      beta1_t = (double *)(&beta1);
      emlrtPopRtStackR2012b(&eb_emlrtRSI, emlrtRootTLSGlobal);
      emlrtPushRtStackR2012b(&fb_emlrtRSI, emlrtRootTLSGlobal);
      Cic0_t = (double *)(&c_y->data[0]);
      emlrtPopRtStackR2012b(&fb_emlrtRSI, emlrtRootTLSGlobal);
      emlrtPushRtStackR2012b(&gb_emlrtRSI, emlrtRootTLSGlobal);
      dgemm(&TRANSA, &TRANSB, &m_t, &n_t, &k_t, alpha1_t, Aia0_t, &lda_t, Bib0_t,
            &ldb_t, beta1_t, Cic0_t, &ldc_t);
      emlrtPopRtStackR2012b(&gb_emlrtRSI, emlrtRootTLSGlobal);
      emlrtPopRtStackR2012b(&o_emlrtRSI, emlrtRootTLSGlobal);
    }

    emlrtPopRtStackR2012b(&m_emlrtRSI, emlrtRootTLSGlobal);
    emlrtPopRtStackR2012b(&i_emlrtRSI, emlrtRootTLSGlobal);
  }

  emxFree_real_T(&a);
  b_emxInit_real_T(&b_x, 1, &e_emlrtRTEI, TRUE);
  i2 = b_x->size[0];
  b_x->size[0] = c_y->size[0];
  emxEnsureCapacity((emxArray__common *)b_x, i2, (int32_T)sizeof(real_T),
                    &e_emlrtRTEI);
  i = c_y->size[0];
  for (i2 = 0; i2 < i; i2++) {
    b_x->data[i2] = c_y->data[i2];
  }

  for (i = 0; i < c_y->size[0]; i++) {
    b_x->data[i] = muDoubleScalarExp(b_x->data[i]);
  }

  i2 = b_x->size[0];
  emxEnsureCapacity((emxArray__common *)b_x, i2, (int32_T)sizeof(real_T),
                    &e_emlrtRTEI);
  i = b_x->size[0];
  for (i2 = 0; i2 < i; i2++) {
    b_x->data[i2]++;
  }

  i2 = c_y->size[0];
  c_y->size[0] = b_x->size[0];
  emxEnsureCapacity((emxArray__common *)c_y, i2, (int32_T)sizeof(real_T),
                    &e_emlrtRTEI);
  i = b_x->size[0];
  for (i2 = 0; i2 < i; i2++) {
    c_y->data[i2] = b_x->data[i2];
  }

  for (i = 0; i < b_x->size[0]; i++) {
    if (b_x->data[i] < 0.0) {
      emlrtPushRtStackR2012b(&e_emlrtRSI, emlrtRootTLSGlobal);
      eml_error();
      emlrtPopRtStackR2012b(&e_emlrtRSI, emlrtRootTLSGlobal);
    }
  }

  for (i = 0; i < b_x->size[0]; i++) {
    c_y->data[i] = muDoubleScalarLog(c_y->data[i]);
  }

  emxFree_real_T(&b_x);
  overflow = FALSE;
  p = FALSE;
  i = 0;
  do {
    exitg1 = 0;
    if (i < 2) {
      if (i + 1 <= 1) {
        i2 = c_y->size[0];
      } else {
        i2 = 1;
      }

      if (i2 != 0) {
        exitg1 = 1;
      } else {
        i++;
      }
    } else {
      p = TRUE;
      exitg1 = 1;
    }
  } while (exitg1 == 0);

  if (!p) {
  } else {
    overflow = TRUE;
  }

  if (!overflow) {
  } else {
    emlrtPushRtStackR2012b(&ib_emlrtRSI, emlrtRootTLSGlobal);
    d_y = NULL;
    m0 = mxCreateCharArray(2, iv4);
    for (i = 0; i < 30; i++) {
      cv5[i] = cv6[i];
    }

    emlrtInitCharArrayR2013a(emlrtRootTLSGlobal, 30, m0, cv5);
    emlrtAssign(&d_y, m0);
    error(message(d_y, &h_emlrtMCI), &i_emlrtMCI);
    emlrtPopRtStackR2012b(&ib_emlrtRSI, emlrtRootTLSGlobal);
  }

  if ((c_y->size[0] == 1) || (c_y->size[0] != 1)) {
    overflow = TRUE;
  } else {
    overflow = FALSE;
  }

  if (overflow) {
  } else {
    emlrtPushRtStackR2012b(&jb_emlrtRSI, emlrtRootTLSGlobal);
    e_y = NULL;
    m0 = mxCreateCharArray(2, iv5);
    for (i = 0; i < 36; i++) {
      cv7[i] = cv8[i];
    }

    emlrtInitCharArrayR2013a(emlrtRootTLSGlobal, 36, m0, cv7);
    emlrtAssign(&e_y, m0);
    error(message(e_y, &j_emlrtMCI), &k_emlrtMCI);
    emlrtPopRtStackR2012b(&jb_emlrtRSI, emlrtRootTLSGlobal);
  }

  if (c_y->size[0] == 0) {
    alpha1 = 0.0;
  } else {
    alpha1 = c_y->data[0];
    emlrtPushRtStackR2012b(&kb_emlrtRSI, emlrtRootTLSGlobal);
    if (2 > c_y->size[0]) {
      overflow = FALSE;
    } else {
      overflow = (c_y->size[0] > 2147483646);
    }

    if (overflow) {
      emlrtPushRtStackR2012b(&t_emlrtRSI, emlrtRootTLSGlobal);
      check_forloop_overflow_error();
      emlrtPopRtStackR2012b(&t_emlrtRSI, emlrtRootTLSGlobal);
    }

    emlrtPopRtStackR2012b(&kb_emlrtRSI, emlrtRootTLSGlobal);
    for (i = 2; i <= c_y->size[0]; i++) {
      alpha1 += c_y->data[i - 1];
    }
  }

  emxFree_real_T(&c_y);
  emxInit_real_T(&b_a, 2, &e_emlrtRTEI, TRUE);
  i2 = b_a->size[0] * b_a->size[1];
  b_a->size[0] = 1;
  emxEnsureCapacity((emxArray__common *)b_a, i2, (int32_T)sizeof(real_T),
                    &e_emlrtRTEI);
  i = x->size[0];
  i2 = b_a->size[0] * b_a->size[1];
  b_a->size[1] = i;
  emxEnsureCapacity((emxArray__common *)b_a, i2, (int32_T)sizeof(real_T),
                    &e_emlrtRTEI);
  i = x->size[0];
  for (i2 = 0; i2 < i; i2++) {
    b_a->data[i2] = x->data[i2];
  }

  emlrtPushRtStackR2012b(&j_emlrtRSI, emlrtRootTLSGlobal);
  if (!(b_a->size[1] == x->size[0])) {
    if ((b_a->size[1] == 1) || (x->size[0] == 1)) {
      emlrtPushRtStackR2012b(&l_emlrtRSI, emlrtRootTLSGlobal);
      f_y = NULL;
      m0 = mxCreateCharArray(2, iv6);
      for (i = 0; i < 45; i++) {
        cv1[i] = cv2[i];
      }

      emlrtInitCharArrayR2013a(emlrtRootTLSGlobal, 45, m0, cv1);
      emlrtAssign(&f_y, m0);
      error(message(f_y, &b_emlrtMCI), &c_emlrtMCI);
      emlrtPopRtStackR2012b(&l_emlrtRSI, emlrtRootTLSGlobal);
    } else {
      emlrtPushRtStackR2012b(&k_emlrtRSI, emlrtRootTLSGlobal);
      g_y = NULL;
      m0 = mxCreateCharArray(2, iv7);
      for (i = 0; i < 21; i++) {
        cv3[i] = cv4[i];
      }

      emlrtInitCharArrayR2013a(emlrtRootTLSGlobal, 21, m0, cv3);
      emlrtAssign(&g_y, m0);
      error(message(g_y, &d_emlrtMCI), &e_emlrtMCI);
      emlrtPopRtStackR2012b(&k_emlrtRSI, emlrtRootTLSGlobal);
    }
  }

  emlrtPopRtStackR2012b(&j_emlrtRSI, emlrtRootTLSGlobal);
  if ((b_a->size[1] == 1) || (x->size[0] == 1)) {
    beta1 = 0.0;
    for (i2 = 0; i2 < b_a->size[1]; i2++) {
      beta1 += b_a->data[b_a->size[0] * i2] * x->data[i2];
    }
  } else {
    emlrtPushRtStackR2012b(&lb_emlrtRSI, emlrtRootTLSGlobal);
    emlrtPushRtStackR2012b(&mb_emlrtRSI, emlrtRootTLSGlobal);
    emlrtPushRtStackR2012b(&nb_emlrtRSI, emlrtRootTLSGlobal);
    if (b_a->size[1] < 1) {
      beta1 = 0.0;
    } else {
      emlrtPushRtStackR2012b(&pb_emlrtRSI, emlrtRootTLSGlobal);
      emlrtPushRtStackR2012b(&sb_emlrtRSI, emlrtRootTLSGlobal);
      emlrtPushRtStackR2012b(&hb_emlrtRSI, emlrtRootTLSGlobal);
      n_t = (ptrdiff_t)(b_a->size[1]);
      emlrtPopRtStackR2012b(&hb_emlrtRSI, emlrtRootTLSGlobal);
      emlrtPopRtStackR2012b(&sb_emlrtRSI, emlrtRootTLSGlobal);
      emlrtPushRtStackR2012b(&tb_emlrtRSI, emlrtRootTLSGlobal);
      emlrtPushRtStackR2012b(&hb_emlrtRSI, emlrtRootTLSGlobal);
      m_t = (ptrdiff_t)(1);
      emlrtPopRtStackR2012b(&hb_emlrtRSI, emlrtRootTLSGlobal);
      emlrtPopRtStackR2012b(&tb_emlrtRSI, emlrtRootTLSGlobal);
      emlrtPushRtStackR2012b(&ub_emlrtRSI, emlrtRootTLSGlobal);
      emlrtPushRtStackR2012b(&hb_emlrtRSI, emlrtRootTLSGlobal);
      k_t = (ptrdiff_t)(1);
      emlrtPopRtStackR2012b(&hb_emlrtRSI, emlrtRootTLSGlobal);
      emlrtPopRtStackR2012b(&ub_emlrtRSI, emlrtRootTLSGlobal);
      emlrtPushRtStackR2012b(&vb_emlrtRSI, emlrtRootTLSGlobal);
      alpha1_t = (double *)(&b_a->data[0]);
      emlrtPopRtStackR2012b(&vb_emlrtRSI, emlrtRootTLSGlobal);
      emlrtPushRtStackR2012b(&wb_emlrtRSI, emlrtRootTLSGlobal);
      Aia0_t = (double *)(&x->data[0]);
      emlrtPopRtStackR2012b(&wb_emlrtRSI, emlrtRootTLSGlobal);
      beta1 = ddot(&n_t, alpha1_t, &m_t, Aia0_t, &k_t);
      emlrtPopRtStackR2012b(&pb_emlrtRSI, emlrtRootTLSGlobal);
    }

    emlrtPopRtStackR2012b(&nb_emlrtRSI, emlrtRootTLSGlobal);
    emlrtPopRtStackR2012b(&mb_emlrtRSI, emlrtRootTLSGlobal);
    emlrtPopRtStackR2012b(&lb_emlrtRSI, emlrtRootTLSGlobal);
  }

  emxFree_real_T(&b_a);
  f = -C * alpha1 - 0.5 * beta1;
  emlrtPopRtStackR2012b(&h_emlrtRSI, emlrtRootTLSGlobal);
  emlrtHeapReferenceStackLeaveFcnR2012b(emlrtRootTLSGlobal);
  return f;
}

/* End of code generation (logpdf.c) */
