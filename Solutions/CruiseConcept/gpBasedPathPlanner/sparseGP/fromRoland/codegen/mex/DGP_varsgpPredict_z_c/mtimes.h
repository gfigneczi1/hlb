/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * mtimes.h
 *
 * Code generation for function 'mtimes'
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
void b_mtimes(const real_T A[140], const real_T B[49], real_T C[140]);

void c_mtimes(const real_T A[140], const real_T B[140], real_T C[400]);

void mtimes(const real_T A[400], const real_T B[400], real_T C[400]);

/* End of code generation (mtimes.h) */
