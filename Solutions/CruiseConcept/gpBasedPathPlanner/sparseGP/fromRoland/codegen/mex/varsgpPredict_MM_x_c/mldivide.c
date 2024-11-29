/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * mldivide.c
 *
 * Code generation for function 'mldivide'
 *
 */

/* Include files */
#include "mldivide.h"
#include "rt_nonfinite.h"
#include "varsgpPredict_MM_x_c_data.h"
#include "warning.h"
#include "xzgetrf.h"
#include <string.h>

/* Variable Definitions */
static emlrtRSInfo
    l_emlrtRSI =
        {
            20,         /* lineNo */
            "mldivide", /* fcnName */
            "/Applications/MATLAB_R2021b.app/toolbox/eml/lib/matlab/ops/"
            "mldivide.m" /* pathName */
};

static emlrtRSInfo
    m_emlrtRSI =
        {
            42,      /* lineNo */
            "mldiv", /* fcnName */
            "/Applications/MATLAB_R2021b.app/toolbox/eml/lib/matlab/ops/"
            "mldivide.m" /* pathName */
};

static emlrtRSInfo n_emlrtRSI =
    {
        67,        /* lineNo */
        "lusolve", /* fcnName */
        "/Applications/MATLAB_R2021b.app/toolbox/eml/eml/+coder/+internal/"
        "lusolve.m" /* pathName */
};

static emlrtRSInfo o_emlrtRSI =
    {
        112,          /* lineNo */
        "lusolveNxN", /* fcnName */
        "/Applications/MATLAB_R2021b.app/toolbox/eml/eml/+coder/+internal/"
        "lusolve.m" /* pathName */
};

static emlrtRSInfo p_emlrtRSI =
    {
        109,          /* lineNo */
        "lusolveNxN", /* fcnName */
        "/Applications/MATLAB_R2021b.app/toolbox/eml/eml/+coder/+internal/"
        "lusolve.m" /* pathName */
};

static emlrtRSInfo q_emlrtRSI =
    {
        124,          /* lineNo */
        "InvAtimesX", /* fcnName */
        "/Applications/MATLAB_R2021b.app/toolbox/eml/eml/+coder/+internal/"
        "lusolve.m" /* pathName */
};

static emlrtRSInfo r_emlrtRSI = {
    26,        /* lineNo */
    "xgetrfs", /* fcnName */
    "/Applications/MATLAB_R2021b.app/toolbox/eml/eml/+coder/+internal/+lapack/"
    "xgetrfs.m" /* pathName */
};

static emlrtRSInfo s_emlrtRSI = {
    27,        /* lineNo */
    "xgetrfs", /* fcnName */
    "/Applications/MATLAB_R2021b.app/toolbox/eml/eml/+coder/+internal/+lapack/"
    "xgetrfs.m" /* pathName */
};

static emlrtRSInfo mb_emlrtRSI =
    {
        90,              /* lineNo */
        "warn_singular", /* fcnName */
        "/Applications/MATLAB_R2021b.app/toolbox/eml/eml/+coder/+internal/"
        "lusolve.m" /* pathName */
};

/* Function Definitions */
void mldivide(const emlrtStack *sp, const real_T A[49], real_T B[49])
{
  emlrtStack b_st;
  emlrtStack c_st;
  emlrtStack d_st;
  emlrtStack e_st;
  emlrtStack f_st;
  emlrtStack g_st;
  emlrtStack st;
  real_T b_A[49];
  real_T temp;
  int32_T ipiv[7];
  int32_T a;
  int32_T b_i;
  int32_T i;
  int32_T i1;
  int32_T info;
  int32_T j;
  int32_T jBcol;
  int32_T k;
  int32_T kAcol;
  st.prev = sp;
  st.tls = sp->tls;
  st.site = &l_emlrtRSI;
  b_st.prev = &st;
  b_st.tls = st.tls;
  c_st.prev = &b_st;
  c_st.tls = b_st.tls;
  d_st.prev = &c_st;
  d_st.tls = c_st.tls;
  e_st.prev = &d_st;
  e_st.tls = d_st.tls;
  f_st.prev = &e_st;
  f_st.tls = e_st.tls;
  g_st.prev = &f_st;
  g_st.tls = f_st.tls;
  b_st.site = &m_emlrtRSI;
  c_st.site = &n_emlrtRSI;
  d_st.site = &p_emlrtRSI;
  e_st.site = &q_emlrtRSI;
  f_st.site = &r_emlrtRSI;
  memcpy(&b_A[0], &A[0], 49U * sizeof(real_T));
  g_st.site = &t_emlrtRSI;
  xzgetrf(&g_st, b_A, ipiv, &info);
  f_st.site = &s_emlrtRSI;
  for (i = 0; i < 6; i++) {
    b_i = ipiv[i];
    if (b_i != i + 1) {
      for (j = 0; j < 7; j++) {
        jBcol = i + 7 * j;
        temp = B[jBcol];
        i1 = (b_i + 7 * j) - 1;
        B[jBcol] = B[i1];
        B[i1] = temp;
      }
    }
  }
  for (j = 0; j < 7; j++) {
    jBcol = 7 * j;
    for (k = 0; k < 7; k++) {
      kAcol = 7 * k;
      b_i = k + jBcol;
      if (B[b_i] != 0.0) {
        a = k + 2;
        for (i = a; i < 8; i++) {
          i1 = (i + jBcol) - 1;
          B[i1] -= B[b_i] * b_A[(i + kAcol) - 1];
        }
      }
    }
  }
  for (j = 0; j < 7; j++) {
    jBcol = 7 * j;
    for (k = 6; k >= 0; k--) {
      kAcol = 7 * k;
      b_i = k + jBcol;
      if (B[b_i] != 0.0) {
        B[b_i] /= b_A[k + kAcol];
        for (i = 0; i < k; i++) {
          i1 = i + jBcol;
          B[i1] -= B[b_i] * b_A[i + kAcol];
        }
      }
    }
  }
  if (info > 0) {
    d_st.site = &o_emlrtRSI;
    e_st.site = &mb_emlrtRSI;
    warning(&e_st);
  }
}

/* End of code generation (mldivide.c) */
