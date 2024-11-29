/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * _coder_DGP_varsgpPredict_y_c_mex.c
 *
 * Code generation for function '_coder_DGP_varsgpPredict_y_c_mex'
 *
 */

/* Include files */
#include "_coder_DGP_varsgpPredict_y_c_mex.h"
#include "DGP_varsgpPredict_y_c_data.h"
#include "DGP_varsgpPredict_y_c_initialize.h"
#include "DGP_varsgpPredict_y_c_terminate.h"
#include "DGP_varsgpPredict_y_c_types.h"
#include "_coder_DGP_varsgpPredict_y_c_api.h"
#include "rt_nonfinite.h"

/* Function Definitions */
void DGP_varsgpPredict_y_c_mexFunction(DGP_varsgpPredict_y_cStackData *SD,
                                       int32_T nlhs, mxArray *plhs[2],
                                       int32_T nrhs, const mxArray *prhs[4])
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
  if (nrhs != 4) {
    emlrtErrMsgIdAndTxt(&st, "EMLRT:runTime:WrongNumberOfInputs", 5, 12, 4, 4,
                        21, "DGP_varsgpPredict_y_c");
  }
  if (nlhs > 2) {
    emlrtErrMsgIdAndTxt(&st, "EMLRT:runTime:TooManyOutputArguments", 3, 4, 21,
                        "DGP_varsgpPredict_y_c");
  }
  /* Call the function. */
  DGP_varsgpPredict_y_c_api(SD, prhs, nlhs, outputs);
  /* Copy over outputs to the caller. */
  if (nlhs < 1) {
    b_nlhs = 1;
  } else {
    b_nlhs = nlhs;
  }
  emlrtReturnArrays(b_nlhs, &plhs[0], &outputs[0]);
}

void mexFunction(int32_T nlhs, mxArray *plhs[], int32_T nrhs,
                 const mxArray *prhs[])
{
  DGP_varsgpPredict_y_cStackData *c_DGP_varsgpPredict_y_cStackDat = NULL;
  c_DGP_varsgpPredict_y_cStackDat =
      (DGP_varsgpPredict_y_cStackData *)emlrtMxCalloc(
          (size_t)1, (size_t)1U * sizeof(DGP_varsgpPredict_y_cStackData));
  mexAtExit(&DGP_varsgpPredict_y_c_atexit);
  /* Module initialization. */
  DGP_varsgpPredict_y_c_initialize();
  /* Dispatch the entry-point. */
  DGP_varsgpPredict_y_c_mexFunction(c_DGP_varsgpPredict_y_cStackDat, nlhs, plhs,
                                    nrhs, prhs);
  /* Module termination. */
  DGP_varsgpPredict_y_c_terminate();
  emlrtMxFree(c_DGP_varsgpPredict_y_cStackDat);
}

emlrtCTX mexFunctionCreateRootTLS(void)
{
  emlrtCreateRootTLSR2021a(&emlrtRootTLSGlobal, &emlrtContextGlobal, NULL, 1,
                           NULL);
  return emlrtRootTLSGlobal;
}

/* End of code generation (_coder_DGP_varsgpPredict_y_c_mex.c) */
