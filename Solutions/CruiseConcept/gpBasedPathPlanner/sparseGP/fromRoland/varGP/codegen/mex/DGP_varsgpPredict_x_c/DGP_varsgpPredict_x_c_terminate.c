/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * DGP_varsgpPredict_x_c_terminate.c
 *
 * Code generation for function 'DGP_varsgpPredict_x_c_terminate'
 *
 */

/* Include files */
#include "DGP_varsgpPredict_x_c_terminate.h"
#include "DGP_varsgpPredict_x_c_data.h"
#include "_coder_DGP_varsgpPredict_x_c_mex.h"
#include "rt_nonfinite.h"

/* Function Definitions */
void DGP_varsgpPredict_x_c_atexit(void)
{
  emlrtStack st = {
      NULL, /* site */
      NULL, /* tls */
      NULL  /* prev */
  };
  mexFunctionCreateRootTLS();
  st.tls = emlrtRootTLSGlobal;
  emlrtEnterRtStackR2012b(&st);
  emlrtLeaveRtStackR2012b(&st);
  emlrtDestroyRootTLS(&emlrtRootTLSGlobal);
  emlrtExitTimeCleanup(&emlrtContextGlobal);
}

void DGP_varsgpPredict_x_c_terminate(void)
{
  emlrtStack st = {
      NULL, /* site */
      NULL, /* tls */
      NULL  /* prev */
  };
  st.tls = emlrtRootTLSGlobal;
  emlrtLeaveRtStackR2012b(&st);
  emlrtDestroyRootTLS(&emlrtRootTLSGlobal);
}

/* End of code generation (DGP_varsgpPredict_x_c_terminate.c) */
