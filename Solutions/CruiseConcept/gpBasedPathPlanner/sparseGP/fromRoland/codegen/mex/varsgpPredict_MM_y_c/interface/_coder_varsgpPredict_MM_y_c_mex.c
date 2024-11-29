/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * _coder_varsgpPredict_MM_y_c_mex.c
 *
 * Code generation for function '_coder_varsgpPredict_MM_y_c_mex'
 *
 */

/* Include files */
#include "_coder_varsgpPredict_MM_y_c_mex.h"
#include "_coder_varsgpPredict_MM_y_c_api.h"
#include "rt_nonfinite.h"
#include "varsgpPredict_MM_y_c_data.h"
#include "varsgpPredict_MM_y_c_initialize.h"
#include "varsgpPredict_MM_y_c_terminate.h"
#include "varsgpPredict_MM_y_c_types.h"

/* Function Definitions */
void mexFunction(int32_T nlhs, mxArray *plhs[], int32_T nrhs,
                 const mxArray *prhs[])
{
  varsgpPredict_MM_y_cStackData *c_varsgpPredict_MM_y_cStackData = NULL;
  c_varsgpPredict_MM_y_cStackData =
      (varsgpPredict_MM_y_cStackData *)emlrtMxCalloc(
          (size_t)1, (size_t)1U * sizeof(varsgpPredict_MM_y_cStackData));
  mexAtExit(&varsgpPredict_MM_y_c_atexit);
  /* Module initialization. */
  varsgpPredict_MM_y_c_initialize();
  /* Dispatch the entry-point. */
  varsgpPredict_MM_y_c_mexFunction(c_varsgpPredict_MM_y_cStackData, nlhs, plhs,
                                   nrhs, prhs);
  /* Module termination. */
  varsgpPredict_MM_y_c_terminate();
  emlrtMxFree(c_varsgpPredict_MM_y_cStackData);
}

emlrtCTX mexFunctionCreateRootTLS(void)
{
  emlrtCreateRootTLSR2021a(&emlrtRootTLSGlobal, &emlrtContextGlobal, NULL, 1,
                           NULL);
  return emlrtRootTLSGlobal;
}

void varsgpPredict_MM_y_c_mexFunction(varsgpPredict_MM_y_cStackData *SD,
                                      int32_T nlhs, mxArray *plhs[2],
                                      int32_T nrhs, const mxArray *prhs[3])
{
  emlrtStack st = {
      NULL, /* site */
      NULL, /* tls */
      NULL  /* prev */
  };
  const mxArray *outputs[2];
  int32_T b_nlhs;
  st.tls = emlrtRootTLSGlobal;
  /* Check for proper number of arguments. */
  if (nrhs != 3) {
    emlrtErrMsgIdAndTxt(&st, "EMLRT:runTime:WrongNumberOfInputs", 5, 12, 3, 4,
                        20, "varsgpPredict_MM_y_c");
  }
  if (nlhs > 2) {
    emlrtErrMsgIdAndTxt(&st, "EMLRT:runTime:TooManyOutputArguments", 3, 4, 20,
                        "varsgpPredict_MM_y_c");
  }
  /* Call the function. */
  varsgpPredict_MM_y_c_api(SD, prhs, nlhs, outputs);
  /* Copy over outputs to the caller. */
  if (nlhs < 1) {
    b_nlhs = 1;
  } else {
    b_nlhs = nlhs;
  }
  emlrtReturnArrays(b_nlhs, &plhs[0], &outputs[0]);
}

/* End of code generation (_coder_varsgpPredict_MM_y_c_mex.c) */
