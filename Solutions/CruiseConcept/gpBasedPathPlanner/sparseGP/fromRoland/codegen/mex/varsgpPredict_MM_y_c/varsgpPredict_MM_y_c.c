/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * varsgpPredict_MM_y_c.c
 *
 * Code generation for function 'varsgpPredict_MM_y_c'
 *
 */

/* Include files */
#include "varsgpPredict_MM_y_c.h"
#include "det.h"
#include "mldivide.h"
#include "mtimes.h"
#include "rt_nonfinite.h"
#include "sumMatrixIncludeNaN.h"
#include "varsgpPredict_MM_y_c_types.h"
#include "blas.h"
#include "mwmathutil.h"
#include <stddef.h>
#include <string.h>

/* Variable Definitions */
static emlrtRSInfo emlrtRSI = {
    50,                     /* lineNo */
    "varsgpPredict_MM_y_c", /* fcnName */
    "/Users/yuhan/Desktop/Dr_Liu_Works/0-文章/0-进行中/2022-Automatica/"
    "Response to the reviewers_Automaica_2023/code for Automatica 2022-v"
    "arGP-new/varGP/varsgpPredict_MM_y_c.m" /* pathName */
};

static emlrtRSInfo b_emlrtRSI = {
    62,                     /* lineNo */
    "varsgpPredict_MM_y_c", /* fcnName */
    "/Users/yuhan/Desktop/Dr_Liu_Works/0-文章/0-进行中/2022-Automatica/"
    "Response to the reviewers_Automaica_2023/code for Automatica 2022-v"
    "arGP-new/varGP/varsgpPredict_MM_y_c.m" /* pathName */
};

static emlrtRSInfo c_emlrtRSI = {
    63,                     /* lineNo */
    "varsgpPredict_MM_y_c", /* fcnName */
    "/Users/yuhan/Desktop/Dr_Liu_Works/0-文章/0-进行中/2022-Automatica/"
    "Response to the reviewers_Automaica_2023/code for Automatica 2022-v"
    "arGP-new/varGP/varsgpPredict_MM_y_c.m" /* pathName */
};

static emlrtRSInfo d_emlrtRSI = {
    64,                     /* lineNo */
    "varsgpPredict_MM_y_c", /* fcnName */
    "/Users/yuhan/Desktop/Dr_Liu_Works/0-文章/0-进行中/2022-Automatica/"
    "Response to the reviewers_Automaica_2023/code for Automatica 2022-v"
    "arGP-new/varGP/varsgpPredict_MM_y_c.m" /* pathName */
};

static emlrtRSInfo e_emlrtRSI = {
    67,                     /* lineNo */
    "varsgpPredict_MM_y_c", /* fcnName */
    "/Users/yuhan/Desktop/Dr_Liu_Works/0-文章/0-进行中/2022-Automatica/"
    "Response to the reviewers_Automaica_2023/code for Automatica 2022-v"
    "arGP-new/varGP/varsgpPredict_MM_y_c.m" /* pathName */
};

static emlrtRSInfo f_emlrtRSI = {
    75,                     /* lineNo */
    "varsgpPredict_MM_y_c", /* fcnName */
    "/Users/yuhan/Desktop/Dr_Liu_Works/0-文章/0-进行中/2022-Automatica/"
    "Response to the reviewers_Automaica_2023/code for Automatica 2022-v"
    "arGP-new/varGP/varsgpPredict_MM_y_c.m" /* pathName */
};

static emlrtRSInfo g_emlrtRSI = {
    77,                     /* lineNo */
    "varsgpPredict_MM_y_c", /* fcnName */
    "/Users/yuhan/Desktop/Dr_Liu_Works/0-文章/0-进行中/2022-Automatica/"
    "Response to the reviewers_Automaica_2023/code for Automatica 2022-v"
    "arGP-new/varGP/varsgpPredict_MM_y_c.m" /* pathName */
};

static emlrtRSInfo h_emlrtRSI = {
    83,                     /* lineNo */
    "varsgpPredict_MM_y_c", /* fcnName */
    "/Users/yuhan/Desktop/Dr_Liu_Works/0-文章/0-进行中/2022-Automatica/"
    "Response to the reviewers_Automaica_2023/code for Automatica 2022-v"
    "arGP-new/varGP/varsgpPredict_MM_y_c.m" /* pathName */
};

static emlrtRTEInfo emlrtRTEI = {
    13,     /* lineNo */
    9,      /* colNo */
    "sqrt", /* fName */
    "/Applications/MATLAB_R2021b.app/toolbox/eml/lib/matlab/elfun/sqrt.m" /* pName
                                                                           */
};

/* Function Definitions */
void varsgpPredict_MM_y_c(const emlrtStack *sp, const struct0_T *model,
                          const real_T m_in[7], const real_T S_in[49],
                          real_T *mu_MM, real_T *var_MM)
{
  ptrdiff_t k_t;
  ptrdiff_t lda_t;
  ptrdiff_t ldb_t;
  ptrdiff_t ldc_t;
  ptrdiff_t m_t;
  ptrdiff_t n_t;
  emlrtStack b_st;
  emlrtStack st;
  real_T A[400];
  real_T b_c[400];
  real_T c[400];
  real_T y[400];
  real_T T[140];
  real_T Zeta_i[140];
  real_T c_y[140];
  real_T d_y[140];
  real_T inp[140];
  real_T R[49];
  real_T c_I[49];
  real_T iLambda[49];
  real_T y_tmp[49];
  real_T beta[20];
  real_T e_y[20];
  real_T log_k[20];
  real_T b_y[7];
  real_T alpha1;
  real_T beta1;
  real_T sigma2f;
  real_T x;
  int32_T b_k;
  int32_T i;
  int32_T k;
  int32_T xoffset;
  char_T TRANSA1;
  char_T TRANSB1;
  int8_T b_I[49];
  int8_T iR_tmp[49];
  st.prev = sp;
  st.tls = sp->tls;
  b_st.prev = &st;
  b_st.tls = st.tls;
  /*  modified from pilco gp1.m  */
  /*  for single dimension */
  sigma2f = muDoubleScalarExp(2.0 * model->GP.logtheta[7]);
  /*  Kmm = Kmm+model.jitter*eye(model.m); */
  /*  L   =  chol(Kmm)';  % K = LL'; */
  /*  invL = inv(L); */
  /*  invKm = invL'*invL; */
  /*  V = L\Kmn; */
  /*   */
  /*  % D = sigma2; */
  /*  V = V./repmat(sqrt(sigma2)',model.m,1); */
  /*  y = model.y./sqrt(sigma2); */
  /*   */
  /*  M = eye(model.m) + V*V'; % M = I+VD-1V' */
  /*  Lm = chol(M)';  % M = LmLm' */
  /*  make predictions  */
  /*  Kmstar = kernel(model.GP, model.Xu, m_in); */
  /*  lst = L\Kmstar; */
  /*  MinvL = Lm*Lm'/L;  */
  /*  LmL = Lm\L'; */
  /*  post_m = MinvL\(V*y); */
  /*  post_S = LmL'*LmL; */
  /*  Kstar = kernel(model.GP, m_in, [], 1);  */
  /*  mustar   = Kmstar'*invKm*post_m; */
  /*  varstar  = Kstar - sum(lst.^2,1)' +
   * diag(Kmstar'*invKm*post_S*invKm*Kmstar); */
  st.site = &emlrtRSI;
  mtimes(model->invKm, model->postS, A);
  st.site = &emlrtRSI;
  mtimes(A, model->invKm, y);
  /*  varstar  = Kstar - diag(Kmstar'*B*Kmstar); */
  /*  MM  */
  for (k = 0; k < 7; k++) {
    for (b_k = 0; b_k < 20; b_k++) {
      xoffset = b_k + 20 * k;
      inp[xoffset] = model->Xu[xoffset] - m_in[k];
    }
  }
  for (b_k = 0; b_k < 20; b_k++) {
    x = 0.0;
    for (i = 0; i < 20; i++) {
      x += model->invKm[b_k + 20 * i] * model->postm[i];
    }
    beta[b_k] = x;
  }
  for (k = 0; k < 7; k++) {
    b_y[k] = muDoubleScalarExp(-2.0 * model->GP.logtheta[k]);
  }
  memset(&iLambda[0], 0, 49U * sizeof(real_T));
  /*  DxD, inversed squared length scale */
  for (xoffset = 0; xoffset < 7; xoffset++) {
    iLambda[xoffset + 7 * xoffset] = b_y[xoffset];
    b_y[xoffset] = muDoubleScalarExp(2.0 * model->GP.logtheta[xoffset]);
  }
  memset(&R[0], 0, 49U * sizeof(real_T));
  for (xoffset = 0; xoffset < 7; xoffset++) {
    R[xoffset + 7 * xoffset] = b_y[xoffset];
  }
  for (b_k = 0; b_k < 49; b_k++) {
    iR_tmp[b_k] = 0;
  }
  for (k = 0; k < 7; k++) {
    iR_tmp[k + 7 * k] = 1;
    for (b_k = 0; b_k < 7; b_k++) {
      x = 0.0;
      for (i = 0; i < 7; i++) {
        x += S_in[k + 7 * i] * iLambda[i + 7 * b_k];
      }
      y_tmp[k + 7 * b_k] = x;
    }
  }
  for (b_k = 0; b_k < 49; b_k++) {
    b_I[b_k] = 0;
  }
  for (k = 0; k < 7; k++) {
    b_I[k + 7 * k] = 1;
  }
  for (b_k = 0; b_k < 49; b_k++) {
    c_I[b_k] = (real_T)b_I[b_k] + y_tmp[b_k];
  }
  st.site = &b_emlrtRSI;
  mldivide(&st, c_I, y_tmp);
  for (b_k = 0; b_k < 49; b_k++) {
    y_tmp[b_k] = (real_T)iR_tmp[b_k] - y_tmp[b_k];
  }
  st.site = &c_emlrtRSI;
  for (b_k = 0; b_k < 7; b_k++) {
    for (i = 0; i < 7; i++) {
      x = 0.0;
      for (xoffset = 0; xoffset < 7; xoffset++) {
        x += iLambda[b_k + 7 * xoffset] * y_tmp[xoffset + 7 * i];
      }
      c_I[b_k + 7 * i] = x;
    }
  }
  b_mtimes(inp, c_I, T);
  st.site = &d_emlrtRSI;
  for (b_k = 0; b_k < 49; b_k++) {
    R[b_k] += S_in[b_k];
  }
  b_st.site = &d_emlrtRSI;
  alpha1 = det(&b_st, R);
  if (alpha1 < 0.0) {
    emlrtErrorWithMessageIdR2018a(
        &st, &emlrtRTEI, "Coder:toolbox:ElFunDomainError",
        "Coder:toolbox:ElFunDomainError", 3, 4, 4, "sqrt");
  }
  alpha1 = muDoubleScalarSqrt(alpha1);
  st.site = &d_emlrtRSI;
  x = b_sumColumnB(*(real_T(*)[7]) & model->GP.logtheta[0]);
  x = muDoubleScalarExp(x);
  x *= sigma2f / alpha1;
  for (b_k = 0; b_k < 140; b_k++) {
    T[b_k] *= inp[b_k];
  }
  memcpy(&log_k[0], &T[0], 20U * sizeof(real_T));
  for (k = 0; k < 6; k++) {
    xoffset = (k + 1) * 20;
    for (b_k = 0; b_k < 20; b_k++) {
      log_k[b_k] += T[xoffset + b_k];
    }
  }
  for (k = 0; k < 20; k++) {
    log_k[k] = muDoubleScalarExp(-log_k[k] / 2.0);
  }
  st.site = &e_emlrtRSI;
  for (b_k = 0; b_k < 20; b_k++) {
    log_k[b_k] = x * log_k[b_k] * beta[b_k];
  }
  *mu_MM = c_sumColumnB(log_k);
  for (k = 0; k < 7; k++) {
    b_y[k] = muDoubleScalarExp(model->GP.logtheta[k]);
    for (b_k = 0; b_k < 20; b_k++) {
      xoffset = b_k + 20 * k;
      T[xoffset] = inp[xoffset] / b_y[k];
    }
  }
  /*  xu-mu/sqrt(iLambda) */
  for (b_k = 0; b_k < 140; b_k++) {
    x = T[b_k];
    x *= x;
    T[b_k] = x;
  }
  memcpy(&log_k[0], &T[0], 20U * sizeof(real_T));
  for (k = 0; k < 6; k++) {
    xoffset = (k + 1) * 20;
    for (b_k = 0; b_k < 20; b_k++) {
      log_k[b_k] += T[xoffset + b_k];
    }
  }
  x = 2.0 * model->GP.logtheta[7];
  for (b_k = 0; b_k < 20; b_k++) {
    log_k[b_k] = x - log_k[b_k] / 2.0;
  }
  /*  DxD */
  for (k = 0; k < 7; k++) {
    x = model->GP.logtheta[k];
    b_y[k] = muDoubleScalarExp(2.0 * x);
    for (b_k = 0; b_k < 20; b_k++) {
      xoffset = b_k + 20 * k;
      Zeta_i[xoffset] = inp[xoffset] / b_y[k];
    }
    x = muDoubleScalarExp(-2.0 * x);
    b_y[k] = x + x;
  }
  memset(&y_tmp[0], 0, 49U * sizeof(real_T));
  for (xoffset = 0; xoffset < 7; xoffset++) {
    y_tmp[xoffset + 7 * xoffset] = b_y[xoffset];
  }
  for (b_k = 0; b_k < 7; b_k++) {
    for (i = 0; i < 7; i++) {
      x = 0.0;
      for (xoffset = 0; xoffset < 7; xoffset++) {
        x += S_in[b_k + 7 * xoffset] * y_tmp[xoffset + 7 * i];
      }
      xoffset = b_k + 7 * i;
      R[xoffset] = x + (real_T)iR_tmp[xoffset];
    }
  }
  st.site = &f_emlrtRSI;
  b_st.site = &f_emlrtRSI;
  alpha1 = det(&b_st, R);
  if (alpha1 < 0.0) {
    emlrtErrorWithMessageIdR2018a(
        &st, &emlrtRTEI, "Coder:toolbox:ElFunDomainError",
        "Coder:toolbox:ElFunDomainError", 3, 4, 4, "sqrt");
  }
  alpha1 = muDoubleScalarSqrt(alpha1);
  x = 1.0 / alpha1;
  st.site = &g_emlrtRSI;
  for (b_k = 0; b_k < 140; b_k++) {
    T[b_k] = -Zeta_i[b_k];
  }
  memcpy(&iLambda[0], &S_in[0], 49U * sizeof(real_T));
  b_st.site = &g_emlrtRSI;
  mldivide(&b_st, R, iLambda);
  for (b_k = 0; b_k < 49; b_k++) {
    iLambda[b_k] /= 2.0;
  }
  /*  maha.m */
  /*  *Summary:* Point-wise squared Mahalanobis distance (a-b)*Q*(a-b)'. */
  /*  Vectors are row-vectors */
  /*  */
  /*     function K = maha(a, b, Q)                          */
  /*  */
  /*  *Input arguments:* */
  /*    */
  /*    a   matrix containing n row vectors                                 [n x
   * D] */
  /*    b   matrix containing n row vectors                                 [n x
   * D] */
  /*    Q   weight matrix. Default: eye(D)                                  [D x
   * D] */
  /*  */
  /*  *Output arguments:* */
  /*   K    point-wise squared distances                                    [n x
   * n] */
  /*  */
  /*  Copyright (C) 2008-2013 by  */
  /*  Marc Deisenroth, Andrew McHutchon, Joe Hall, and Carl Edward Rasmussen. */
  /*  */
  /*  Last modified: 2013-03-21 */
  /*  Code */
  b_mtimes(Zeta_i, iLambda, inp);
  b_mtimes(T, iLambda, c_y);
  for (b_k = 0; b_k < 140; b_k++) {
    d_y[b_k] = 2.0 * inp[b_k];
  }
  TRANSB1 = 'T';
  TRANSA1 = 'N';
  alpha1 = 1.0;
  beta1 = 0.0;
  m_t = (ptrdiff_t)20;
  n_t = (ptrdiff_t)20;
  k_t = (ptrdiff_t)7;
  lda_t = (ptrdiff_t)20;
  ldb_t = (ptrdiff_t)20;
  ldc_t = (ptrdiff_t)20;
  dgemm(&TRANSA1, &TRANSB1, &m_t, &n_t, &k_t, &alpha1, &d_y[0], &lda_t, &T[0],
        &ldb_t, &beta1, &A[0], &ldc_t);
  for (k = 0; k < 20; k++) {
    for (b_k = 0; b_k < 20; b_k++) {
      c[b_k + 20 * k] = log_k[b_k] + log_k[k];
    }
  }
  for (b_k = 0; b_k < 140; b_k++) {
    inp[b_k] *= Zeta_i[b_k];
  }
  memcpy(&log_k[0], &inp[0], 20U * sizeof(real_T));
  for (k = 0; k < 6; k++) {
    xoffset = (k + 1) * 20;
    for (b_k = 0; b_k < 20; b_k++) {
      log_k[b_k] += inp[xoffset + b_k];
    }
  }
  for (b_k = 0; b_k < 140; b_k++) {
    c_y[b_k] *= T[b_k];
  }
  memcpy(&e_y[0], &c_y[0], 20U * sizeof(real_T));
  for (k = 0; k < 6; k++) {
    xoffset = (k + 1) * 20;
    for (b_k = 0; b_k < 20; b_k++) {
      e_y[b_k] += c_y[xoffset + b_k];
    }
  }
  for (k = 0; k < 20; k++) {
    for (b_k = 0; b_k < 20; b_k++) {
      b_c[b_k + 20 * k] = log_k[b_k] + e_y[k];
    }
  }
  for (k = 0; k < 400; k++) {
    c[k] = muDoubleScalarExp(c[k] + (b_c[k] - A[k]));
  }
  for (b_k = 0; b_k < 20; b_k++) {
    for (i = 0; i < 20; i++) {
      xoffset = i + 20 * b_k;
      A[xoffset] = beta[i] * beta[b_k] - (model->invKm[xoffset] - y[xoffset]);
    }
  }
  for (b_k = 0; b_k < 400; b_k++) {
    A[b_k] *= x * c[b_k];
  }
  st.site = &h_emlrtRSI;
  for (xoffset = 0; xoffset < 20; xoffset++) {
    log_k[xoffset] = sumColumnB(A, xoffset + 1);
  }
  st.site = &h_emlrtRSI;
  *var_MM = c_sumColumnB(log_k);
  *var_MM += sigma2f;
  *var_MM -= *mu_MM * *mu_MM;
  /*  if var_MM > sigma2f */
  /*      var_MM = sigma2f; */
  /*  end */
}

/* End of code generation (varsgpPredict_MM_y_c.c) */
