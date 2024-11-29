/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * _coder_varsgpPredict_MM_z_c_mex.h
 *
 * Code generation for function '_coder_varsgpPredict_MM_z_c_mex'
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
MEXFUNCTION_LINKAGE void mexFunction(int32_T nlhs, mxArray *plhs[],
                                     int32_T nrhs, const mxArray *prhs[]);

emlrtCTX mexFunctionCreateRootTLS(void);

void varsgpPredict_MM_z_c_mexFunction(varsgpPredict_MM_z_cStackData *SD,
                                      int32_T nlhs, mxArray *plhs[2],
                                      int32_T nrhs, const mxArray *prhs[3]);

/* End of code generation (_coder_varsgpPredict_MM_z_c_mex.h) */
