/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * xzgetrf.c
 *
 * Code generation for function 'xzgetrf'
 *
 */

/* Include files */
#include "xzgetrf.h"
#include "eml_int_forloop_overflow_check.h"
#include "rt_nonfinite.h"
#include "mwmathutil.h"

/* Variable Definitions */
static emlrtRSInfo fb_emlrtRSI = {
    36,        /* lineNo */
    "xzgetrf", /* fcnName */
    "/Applications/MATLAB_R2021b.app/toolbox/eml/eml/+coder/+internal/"
    "+reflapack/xzgetrf.m" /* pathName */
};

static emlrtRSInfo gb_emlrtRSI = {
    50,        /* lineNo */
    "xzgetrf", /* fcnName */
    "/Applications/MATLAB_R2021b.app/toolbox/eml/eml/+coder/+internal/"
    "+reflapack/xzgetrf.m" /* pathName */
};

static emlrtRSInfo hb_emlrtRSI = {
    58,        /* lineNo */
    "xzgetrf", /* fcnName */
    "/Applications/MATLAB_R2021b.app/toolbox/eml/eml/+coder/+internal/"
    "+reflapack/xzgetrf.m" /* pathName */
};

static emlrtRSInfo ib_emlrtRSI = {
    23,       /* lineNo */
    "ixamax", /* fcnName */
    "/Applications/MATLAB_R2021b.app/toolbox/eml/eml/+coder/+internal/+blas/"
    "ixamax.m" /* pathName */
};

static emlrtRSInfo jb_emlrtRSI = {
    24,       /* lineNo */
    "ixamax", /* fcnName */
    "/Applications/MATLAB_R2021b.app/toolbox/eml/eml/+coder/+internal/+refblas/"
    "ixamax.m" /* pathName */
};

static emlrtRSInfo kb_emlrtRSI = {
    21,                               /* lineNo */
    "eml_int_forloop_overflow_check", /* fcnName */
    "/Applications/MATLAB_R2021b.app/toolbox/eml/lib/matlab/eml/"
    "eml_int_forloop_overflow_check.m" /* pathName */
};

static emlrtRSInfo lb_emlrtRSI = {
    45,      /* lineNo */
    "xgeru", /* fcnName */
    "/Applications/MATLAB_R2021b.app/toolbox/eml/eml/+coder/+internal/+blas/"
    "xgeru.m" /* pathName */
};

static emlrtRSInfo mb_emlrtRSI = {
    45,     /* lineNo */
    "xger", /* fcnName */
    "/Applications/MATLAB_R2021b.app/toolbox/eml/eml/+coder/+internal/+blas/"
    "xger.m" /* pathName */
};

static emlrtRSInfo nb_emlrtRSI = {
    15,     /* lineNo */
    "xger", /* fcnName */
    "/Applications/MATLAB_R2021b.app/toolbox/eml/eml/+coder/+internal/+refblas/"
    "xger.m" /* pathName */
};

static emlrtRSInfo ob_emlrtRSI = {
    41,      /* lineNo */
    "xgerx", /* fcnName */
    "/Applications/MATLAB_R2021b.app/toolbox/eml/eml/+coder/+internal/+refblas/"
    "xgerx.m" /* pathName */
};

static emlrtRSInfo pb_emlrtRSI = {
    54,      /* lineNo */
    "xgerx", /* fcnName */
    "/Applications/MATLAB_R2021b.app/toolbox/eml/eml/+coder/+internal/+refblas/"
    "xgerx.m" /* pathName */
};

/* Function Definitions */
void xzgetrf(const emlrtStack *sp, real_T A[49], int32_T ipiv[7], int32_T *info)
{
  emlrtStack b_st;
  emlrtStack c_st;
  emlrtStack d_st;
  emlrtStack e_st;
  emlrtStack f_st;
  emlrtStack st;
  real_T s;
  real_T smax;
  int32_T b_tmp;
  int32_T j;
  int32_T jA;
  int32_T jp1j;
  int32_T k;
  int32_T n;
  int32_T temp_tmp;
  st.prev = sp;
  st.tls = sp->tls;
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
  for (n = 0; n < 7; n++) {
    ipiv[n] = n + 1;
  }
  *info = 0;
  for (j = 0; j < 6; j++) {
    b_tmp = j << 3;
    jp1j = b_tmp + 2;
    n = 7 - j;
    st.site = &fb_emlrtRSI;
    b_st.site = &ib_emlrtRSI;
    jA = 0;
    smax = muDoubleScalarAbs(A[b_tmp]);
    c_st.site = &jb_emlrtRSI;
    for (k = 2; k <= n; k++) {
      s = muDoubleScalarAbs(A[(b_tmp + k) - 1]);
      if (s > smax) {
        jA = k - 1;
        smax = s;
      }
    }
    if (A[b_tmp + jA] != 0.0) {
      if (jA != 0) {
        jA += j;
        ipiv[j] = jA + 1;
        for (k = 0; k < 7; k++) {
          temp_tmp = j + k * 7;
          smax = A[temp_tmp];
          n = jA + k * 7;
          A[temp_tmp] = A[n];
          A[n] = smax;
        }
      }
      temp_tmp = (b_tmp - j) + 7;
      st.site = &gb_emlrtRSI;
      for (jA = jp1j; jA <= temp_tmp; jA++) {
        A[jA - 1] /= A[b_tmp];
      }
    } else {
      *info = j + 1;
    }
    n = 5 - j;
    st.site = &hb_emlrtRSI;
    b_st.site = &lb_emlrtRSI;
    c_st.site = &mb_emlrtRSI;
    d_st.site = &nb_emlrtRSI;
    jA = b_tmp + 9;
    e_st.site = &ob_emlrtRSI;
    for (jp1j = 0; jp1j <= n; jp1j++) {
      smax = A[(b_tmp + jp1j * 7) + 7];
      if (smax != 0.0) {
        temp_tmp = (jA - j) + 5;
        e_st.site = &pb_emlrtRSI;
        if ((jA <= temp_tmp) && (temp_tmp > 2147483646)) {
          f_st.site = &kb_emlrtRSI;
          check_forloop_overflow_error(&f_st);
        }
        for (k = jA; k <= temp_tmp; k++) {
          A[k - 1] += A[((b_tmp + k) - jA) + 1] * -smax;
        }
      }
      jA += 7;
    }
  }
  if ((*info == 0) && (!(A[48] != 0.0))) {
    *info = 7;
  }
}

/* End of code generation (xzgetrf.c) */
