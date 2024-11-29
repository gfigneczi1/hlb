/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * DGP_varsgpPredict_z_c_initialize.c
 *
 * Code generation for function 'DGP_varsgpPredict_z_c_initialize'
 *
 */

/* Include files */
#include "DGP_varsgpPredict_z_c_initialize.h"
#include "DGP_varsgpPredict_z_c_data.h"
#include "_coder_DGP_varsgpPredict_z_c_mex.h"
#include "rt_nonfinite.h"

/* Function Definitions */
void DGP_varsgpPredict_z_c_initialize(void)
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

/* End of code generation (DGP_varsgpPredict_z_c_initialize.c) */
