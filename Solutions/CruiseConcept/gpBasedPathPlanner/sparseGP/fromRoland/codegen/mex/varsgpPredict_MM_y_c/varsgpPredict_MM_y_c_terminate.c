/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * varsgpPredict_MM_y_c_terminate.c
 *
 * Code generation for function 'varsgpPredict_MM_y_c_terminate'
 *
 */

/* Include files */
#include "varsgpPredict_MM_y_c_terminate.h"
#include "_coder_varsgpPredict_MM_y_c_mex.h"
#include "rt_nonfinite.h"
#include "varsgpPredict_MM_y_c_data.h"

/* Function Definitions */
void varsgpPredict_MM_y_c_atexit(void)
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

void varsgpPredict_MM_y_c_terminate(void)
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

/* End of code generation (varsgpPredict_MM_y_c_terminate.c) */
