/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * xzgetrf.h
 *
 * Code generation for function 'xzgetrf'
 *
 */

#pragma once

/* Include files */
#include "rtwtypes.h"
#include "emlrt.h"
#include "mex.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* Function Declarations */
void xzgetrf(const emlrtStack *sp, real_T A[49], int32_T ipiv[7],
             int32_T *info);

/* End of code generation (xzgetrf.h) */
