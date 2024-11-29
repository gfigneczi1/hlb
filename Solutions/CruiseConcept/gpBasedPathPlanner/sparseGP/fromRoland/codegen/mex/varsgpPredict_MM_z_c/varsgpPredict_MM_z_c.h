/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * varsgpPredict_MM_z_c.h
 *
 * Code generation for function 'varsgpPredict_MM_z_c'
 *
 */

#pragma once

/* Include files */
#include "rtwtypes.h"
#include "varsgpPredict_MM_z_c_types.h"
#include "emlrt.h"
#include "mex.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* Function Declarations */
void varsgpPredict_MM_z_c(const emlrtStack *sp, const struct0_T *model,
                          const real_T m_in[7], const real_T S_in[49],
                          real_T *mu_MM, real_T *var_MM);

/* End of code generation (varsgpPredict_MM_z_c.h) */
