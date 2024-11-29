/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * det.c
 *
 * Code generation for function 'det'
 *
 */

/* Include files */
#include "det.h"
#include "rt_nonfinite.h"
#include "varsgpPredict_MM_z_c_data.h"
#include "xzgetrf.h"
#include <string.h>

/* Variable Definitions */
static emlrtRSInfo nb_emlrtRSI = {
    21,    /* lineNo */
    "det", /* fcnName */
    "/Applications/MATLAB_R2021b.app/toolbox/eml/lib/matlab/matfun/det.m" /* pathName
                                                                           */
};

/* Function Definitions */
real_T det(const emlrtStack *sp, const real_T x[49])
{
  emlrtStack b_st;
  emlrtStack st;
  real_T b_x[49];
  real_T y;
  int32_T ipiv[7];
  int32_T info;
  boolean_T isodd;
  st.prev = sp;
  st.tls = sp->tls;
  st.site = &nb_emlrtRSI;
  b_st.prev = &st;
  b_st.tls = st.tls;
  memcpy(&b_x[0], &x[0], 49U * sizeof(real_T));
  b_st.site = &t_emlrtRSI;
  xzgetrf(&b_st, b_x, ipiv, &info);
  y = b_x[0];
  isodd = false;
  for (info = 0; info < 6; info++) {
    y *= b_x[(info + 7 * (info + 1)) + 1];
    if (ipiv[info] > info + 1) {
      isodd = !isodd;
    }
  }
  if (isodd) {
    y = -y;
  }
  return y;
}

/* End of code generation (det.c) */
