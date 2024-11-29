/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * varsgpPredict_MM_y_c_initialize.c
 *
 * Code generation for function 'varsgpPredict_MM_y_c_initialize'
 *
 */

/* Include files */
#include "varsgpPredict_MM_y_c_initialize.h"
#include "_coder_varsgpPredict_MM_y_c_mex.h"
#include "rt_nonfinite.h"
#include "varsgpPredict_MM_y_c_data.h"

/* Function Definitions */
void varsgpPredict_MM_y_c_initialize(void)
{
  static const volatile char_T *emlrtBreakCheckR2012bFlagVar = NULL;
  emlrtStack st = {
      NULL, /* site */
      NULL, /* tls */
      NULL  /* prev */
  };
  mex_InitInfAndNan();
  mexFunctionCreateRootTLS();
  emlrtBreakCheckR2012bFlagVar = emlrtGetBreakCheckFlagAddressR2012b();
  st.tls = emlrtRootTLSGlobal;
  emlrtClearAllocCountR2012b(&st, false, 0U, NULL);
  emlrtEnterRtStackR2012b(&st);
  emlrtFirstTimeR2012b(emlrtRootTLSGlobal);
}

/* End of code generation (varsgpPredict_MM_y_c_initialize.c) */
