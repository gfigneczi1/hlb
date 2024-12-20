/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * sumMatrixIncludeNaN.c
 *
 * Code generation for function 'sumMatrixIncludeNaN'
 *
 */

/* Include files */
#include "sumMatrixIncludeNaN.h"
#include "rt_nonfinite.h"

/* Function Definitions */
real_T b_sumColumnB(const real_T x[7])
{
  real_T y;
  int32_T k;
  y = x[0];
  for (k = 0; k < 6; k++) {
    y += x[k + 1];
  }
  return y;
}

real_T c_sumColumnB(const real_T x[20])
{
  real_T y;
  int32_T k;
  y = x[0];
  for (k = 0; k < 19; k++) {
    y += x[k + 1];
  }
  return y;
}

real_T sumColumnB(const real_T x[400], int32_T col)
{
  real_T y;
  int32_T i0;
  int32_T k;
  i0 = (col - 1) * 20;
  y = x[i0];
  for (k = 0; k < 19; k++) {
    y += x[(i0 + k) + 1];
  }
  return y;
}

/* End of code generation (sumMatrixIncludeNaN.c) */
