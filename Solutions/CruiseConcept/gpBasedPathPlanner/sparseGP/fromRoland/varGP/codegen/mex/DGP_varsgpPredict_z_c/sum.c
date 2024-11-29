/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * sum.c
 *
 * Code generation for function 'sum'
 *
 */

/* Include files */
#include "sum.h"
#include "rt_nonfinite.h"
#include <string.h>

/* Function Definitions */
void sum(const real_T x[140], real_T y[20])
{
  int32_T k;
  int32_T xj;
  int32_T xoffset;
  memcpy(&y[0], &x[0], 20U * sizeof(real_T));
  for (k = 0; k < 6; k++) {
    xoffset = (k + 1) * 20;
    for (xj = 0; xj < 20; xj++) {
      y[xj] += x[xoffset + xj];
    }
  }
}

/* End of code generation (sum.c) */
