/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * DGP_varsgpPredict_y_c.c
 *
 * Code generation for function 'DGP_varsgpPredict_y_c'
 *
 */

/* Include files */
#include "DGP_varsgpPredict_y_c.h"
#include "DGP_varsgpPredict_y_c_types.h"
#include "det.h"
#include "mldivide.h"
#include "mtimes.h"
#include "rt_nonfinite.h"
#include "sum.h"
#include "sumMatrixIncludeNaN.h"
#include "mwmathutil.h"
#include <string.h>

/* Variable Definitions */
static emlrtRSInfo emlrtRSI = {
    67,                      /* lineNo */
    "DGP_varsgpPredict_y_c", /* fcnName */
    "/Users/yuhan/Desktop/Dr_Liu_Works/0-文章/0-进行中/2022-Automatica/"
    "Response to the reviewers_Automaica_2023/code for Automatica 2022-v"
    "arGP-new/varGP/DGP_varsgpPredict_y_c.m" /* pathName */
};

static emlrtRSInfo b_emlrtRSI = {
    68,                      /* lineNo */
    "DGP_varsgpPredict_y_c", /* fcnName */
    "/Users/yuhan/Desktop/Dr_Liu_Works/0-文章/0-进行中/2022-Automatica/"
    "Response to the reviewers_Automaica_2023/code for Automatica 2022-v"
    "arGP-new/varGP/DGP_varsgpPredict_y_c.m" /* pathName */
};

static emlrtRSInfo c_emlrtRSI = {
    83,                      /* lineNo */
    "DGP_varsgpPredict_y_c", /* fcnName */
    "/Users/yuhan/Desktop/Dr_Liu_Works/0-文章/0-进行中/2022-Automatica/"
    "Response to the reviewers_Automaica_2023/code for Automatica 2022-v"
    "arGP-new/varGP/DGP_varsgpPredict_y_c.m" /* pathName */
};

static emlrtRSInfo d_emlrtRSI = {
    84,                      /* lineNo */
    "DGP_varsgpPredict_y_c", /* fcnName */
    "/Users/yuhan/Desktop/Dr_Liu_Works/0-文章/0-进行中/2022-Automatica/"
    "Response to the reviewers_Automaica_2023/code for Automatica 2022-v"
    "arGP-new/varGP/DGP_varsgpPredict_y_c.m" /* pathName */
};

static emlrtRSInfo e_emlrtRSI = {
    85,                      /* lineNo */
    "DGP_varsgpPredict_y_c", /* fcnName */
    "/Users/yuhan/Desktop/Dr_Liu_Works/0-文章/0-进行中/2022-Automatica/"
    "Response to the reviewers_Automaica_2023/code for Automatica 2022-v"
    "arGP-new/varGP/DGP_varsgpPredict_y_c.m" /* pathName */
};

static emlrtRSInfo f_emlrtRSI = {
    95,                      /* lineNo */
    "DGP_varsgpPredict_y_c", /* fcnName */
    "/Users/yuhan/Desktop/Dr_Liu_Works/0-文章/0-进行中/2022-Automatica/"
    "Response to the reviewers_Automaica_2023/code for Automatica 2022-v"
    "arGP-new/varGP/DGP_varsgpPredict_y_c.m" /* pathName */
};

static emlrtRSInfo g_emlrtRSI = {
    96,                      /* lineNo */
    "DGP_varsgpPredict_y_c", /* fcnName */
    "/Users/yuhan/Desktop/Dr_Liu_Works/0-文章/0-进行中/2022-Automatica/"
    "Response to the reviewers_Automaica_2023/code for Automatica 2022-v"
    "arGP-new/varGP/DGP_varsgpPredict_y_c.m" /* pathName */
};

static emlrtRSInfo h_emlrtRSI = {
    97,                      /* lineNo */
    "DGP_varsgpPredict_y_c", /* fcnName */
    "/Users/yuhan/Desktop/Dr_Liu_Works/0-文章/0-进行中/2022-Automatica/"
    "Response to the reviewers_Automaica_2023/code for Automatica 2022-v"
    "arGP-new/varGP/DGP_varsgpPredict_y_c.m" /* pathName */
};

static emlrtRSInfo i_emlrtRSI = {
    102,                     /* lineNo */
    "DGP_varsgpPredict_y_c", /* fcnName */
    "/Users/yuhan/Desktop/Dr_Liu_Works/0-文章/0-进行中/2022-Automatica/"
    "Response to the reviewers_Automaica_2023/code for Automatica 2022-v"
    "arGP-new/varGP/DGP_varsgpPredict_y_c.m" /* pathName */
};

static emlrtRSInfo j_emlrtRSI = {
    112,                     /* lineNo */
    "DGP_varsgpPredict_y_c", /* fcnName */
    "/Users/yuhan/Desktop/Dr_Liu_Works/0-文章/0-进行中/2022-Automatica/"
    "Response to the reviewers_Automaica_2023/code for Automatica 2022-v"
    "arGP-new/varGP/DGP_varsgpPredict_y_c.m" /* pathName */
};

static emlrtRSInfo k_emlrtRSI = {
    114,                     /* lineNo */
    "DGP_varsgpPredict_y_c", /* fcnName */
    "/Users/yuhan/Desktop/Dr_Liu_Works/0-文章/0-进行中/2022-Automatica/"
    "Response to the reviewers_Automaica_2023/code for Automatica 2022-v"
    "arGP-new/varGP/DGP_varsgpPredict_y_c.m" /* pathName */
};

static emlrtRSInfo l_emlrtRSI = {
    120,                     /* lineNo */
    "DGP_varsgpPredict_y_c", /* fcnName */
    "/Users/yuhan/Desktop/Dr_Liu_Works/0-文章/0-进行中/2022-Automatica/"
    "Response to the reviewers_Automaica_2023/code for Automatica 2022-v"
    "arGP-new/varGP/DGP_varsgpPredict_y_c.m" /* pathName */
};

static emlrtRSInfo m_emlrtRSI = {
    132,                     /* lineNo */
    "DGP_varsgpPredict_y_c", /* fcnName */
    "/Users/yuhan/Desktop/Dr_Liu_Works/0-文章/0-进行中/2022-Automatica/"
    "Response to the reviewers_Automaica_2023/code for Automatica 2022-v"
    "arGP-new/varGP/DGP_varsgpPredict_y_c.m" /* pathName */
};

static emlrtRSInfo n_emlrtRSI = {
    134,                     /* lineNo */
    "DGP_varsgpPredict_y_c", /* fcnName */
    "/Users/yuhan/Desktop/Dr_Liu_Works/0-文章/0-进行中/2022-Automatica/"
    "Response to the reviewers_Automaica_2023/code for Automatica 2022-v"
    "arGP-new/varGP/DGP_varsgpPredict_y_c.m" /* pathName */
};

static emlrtRSInfo o_emlrtRSI = {
    140,                     /* lineNo */
    "DGP_varsgpPredict_y_c", /* fcnName */
    "/Users/yuhan/Desktop/Dr_Liu_Works/0-文章/0-进行中/2022-Automatica/"
    "Response to the reviewers_Automaica_2023/code for Automatica 2022-v"
    "arGP-new/varGP/DGP_varsgpPredict_y_c.m" /* pathName */
};

static emlrtRSInfo p_emlrtRSI = {
    150,                     /* lineNo */
    "DGP_varsgpPredict_y_c", /* fcnName */
    "/Users/yuhan/Desktop/Dr_Liu_Works/0-文章/0-进行中/2022-Automatica/"
    "Response to the reviewers_Automaica_2023/code for Automatica 2022-v"
    "arGP-new/varGP/DGP_varsgpPredict_y_c.m" /* pathName */
};

static emlrtRSInfo q_emlrtRSI = {
    152,                     /* lineNo */
    "DGP_varsgpPredict_y_c", /* fcnName */
    "/Users/yuhan/Desktop/Dr_Liu_Works/0-文章/0-进行中/2022-Automatica/"
    "Response to the reviewers_Automaica_2023/code for Automatica 2022-v"
    "arGP-new/varGP/DGP_varsgpPredict_y_c.m" /* pathName */
};

static emlrtRSInfo r_emlrtRSI = {
    156,                     /* lineNo */
    "DGP_varsgpPredict_y_c", /* fcnName */
    "/Users/yuhan/Desktop/Dr_Liu_Works/0-文章/0-进行中/2022-Automatica/"
    "Response to the reviewers_Automaica_2023/code for Automatica 2022-v"
    "arGP-new/varGP/DGP_varsgpPredict_y_c.m" /* pathName */
};

static emlrtRTEInfo b_emlrtRTEI = {
    13,     /* lineNo */
    9,      /* colNo */
    "sqrt", /* fName */
    "/Applications/MATLAB_R2021b.app/toolbox/eml/lib/matlab/elfun/sqrt.m" /* pName
                                                                           */
};

/* Function Definitions */
void DGP_varsgpPredict_y_c(const emlrtStack *sp, const struct0_T *model,
                           const struct3_T *Smodel, const real_T m_in[7],
                           const real_T S_in[49], real_T *mu_MM, real_T *var_MM)
{
  emlrtStack b_st;
  emlrtStack st;
  real_T A[400];
  real_T b_c[400];
  real_T b_y[400];
  real_T c[400];
  real_T y[400];
  real_T T[140];
  real_T Zeta_i[140];
  real_T aQ[140];
  real_T dv[140];
  real_T e_y[140];
  real_T inp[140];
  real_T inp_s[140];
  real_T c_I[49];
  real_T d[49];
  real_T iLambda[49];
  real_T iLambda_s[49];
  real_T y_tmp[49];
  real_T b_x[20];
  real_T beta[20];
  real_T beta_s[20];
  real_T log_k[20];
  real_T log_ks[20];
  real_T x[20];
  real_T d_y[7];
  real_T a;
  real_T c_y;
  real_T sigma2f;
  real_T sigma2f_s;
  real_T var_long;
  real_T var_short;
  int32_T b_k;
  int32_T i;
  int32_T inp_tmp;
  int32_T k;
  int8_T b_I[49];
  int8_T iR_tmp[49];
  st.prev = sp;
  st.tls = sp->tls;
  b_st.prev = &st;
  b_st.tls = st.tls;
  /*  modified from pilco gp1.m  */
  /*  for single dimension */
  /*  sigma2 = exp(2*model.Likelihood.logtheta); */
  /*  sigma2_s = exp(2*Smodel.Likelihood.logtheta); */
  sigma2f = muDoubleScalarExp(2.0 * model->GP.logtheta[7]);
  sigma2f_s = muDoubleScalarExp(2.0 * Smodel->GP.logtheta[7]);
  /*  long-term part */
  /*  V = L\Kmn; */
  /*   */
  /*  V = V./repmat(sqrt(sigma2)',model.m,1); */
  /*  y = model.y./sqrt(sigma2); */
  /*   */
  /*  M = eye(model.m) + V*V'; % M = I+VD-1V' */
  /*  Lm = chol(M)';  % M = LmLm' */
  /*  short-term part */
  /*  V_s = L_s\Kmn_s; */
  /*   */
  /*  % D = sigma2; */
  /*  V_s = V_s./repmat(sqrt(sigma2_s)',Smodel.m,1); */
  /*  y_s = Smodel.y./sqrt(sigma2_s); */
  /*   */
  /*  M_s = eye(Smodel.m) + V_s*V_s'; % M = I+VD-1V' */
  /*  Lm_s = chol(M_s)';  % M = LmLm' */
  /*  make predictions  */
  /*  Kmstar = kernel(model.GP, model.Xu, m_in); */
  /*  Kmstar_s = kernel(Smodel.GP, Smodel.Xu, m_in); */
  /*  lst = L\Kmstar; */
  /*  Kstar = kernel(model.GP, m_in, [], 1);  */
  /*  Kstar_s = kernel(Smodel.GP, m_in, [], 1);  */
  /*  mustar   = Kmstar'*invKm*post_m +  Kmstar_s'*invKm_s*Spost_m ; */
  /*  varstar  = Kstar - sum(lst.^2,1)' +
   * diag(Kmstar'*invKm*post_S*invKm*Kmstar); */
  st.site = &emlrtRSI;
  mtimes(model->invKm, model->postS, A);
  st.site = &emlrtRSI;
  mtimes(A, model->invKm, y);
  st.site = &b_emlrtRSI;
  mtimes(Smodel->invKm, Smodel->postS, A);
  st.site = &b_emlrtRSI;
  mtimes(A, Smodel->invKm, b_y);
  /*  long_var = Kstar - diag(Kmstar'*B*Kmstar); */
  /*  short_var = Kstar_s - diag(Kmstar_s'*B_s*Kmstar_s); */
  /*  varstar  = long_var + short_var; */
  /*  MM  */
  /*  long-term mean */
  for (k = 0; k < 7; k++) {
    for (b_k = 0; b_k < 20; b_k++) {
      inp_tmp = b_k + 20 * k;
      inp[inp_tmp] = model->Xu[inp_tmp] - m_in[k];
    }
  }
  for (b_k = 0; b_k < 20; b_k++) {
    c_y = 0.0;
    for (i = 0; i < 20; i++) {
      c_y += model->invKm[b_k + 20 * i] * model->postm[i];
    }
    beta[b_k] = c_y;
  }
  for (k = 0; k < 7; k++) {
    d_y[k] = muDoubleScalarExp(-2.0 * model->GP.logtheta[k]);
  }
  memset(&iLambda[0], 0, 49U * sizeof(real_T));
  /*  DxD, inversed squared length scale */
  for (inp_tmp = 0; inp_tmp < 7; inp_tmp++) {
    iLambda[inp_tmp + 7 * inp_tmp] = d_y[inp_tmp];
    d_y[inp_tmp] = muDoubleScalarExp(2.0 * model->GP.logtheta[inp_tmp]);
  }
  memset(&d[0], 0, 49U * sizeof(real_T));
  for (inp_tmp = 0; inp_tmp < 7; inp_tmp++) {
    d[inp_tmp + 7 * inp_tmp] = d_y[inp_tmp];
  }
  for (b_k = 0; b_k < 49; b_k++) {
    iR_tmp[b_k] = 0;
  }
  for (k = 0; k < 7; k++) {
    iR_tmp[k + 7 * k] = 1;
    for (b_k = 0; b_k < 7; b_k++) {
      c_y = 0.0;
      for (i = 0; i < 7; i++) {
        c_y += S_in[k + 7 * i] * iLambda[i + 7 * b_k];
      }
      y_tmp[k + 7 * b_k] = c_y;
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
  st.site = &c_emlrtRSI;
  mldivide(&st, c_I, y_tmp);
  for (b_k = 0; b_k < 49; b_k++) {
    y_tmp[b_k] = (real_T)iR_tmp[b_k] - y_tmp[b_k];
  }
  st.site = &d_emlrtRSI;
  for (b_k = 0; b_k < 7; b_k++) {
    for (i = 0; i < 7; i++) {
      c_y = 0.0;
      for (inp_tmp = 0; inp_tmp < 7; inp_tmp++) {
        c_y += iLambda[b_k + 7 * inp_tmp] * y_tmp[inp_tmp + 7 * i];
      }
      c_I[b_k + 7 * i] = c_y;
    }
  }
  b_mtimes(inp, c_I, T);
  st.site = &e_emlrtRSI;
  var_short = b_sumColumnB(*(real_T(*)[7]) & model->GP.logtheta[0]);
  var_short = muDoubleScalarExp(var_short);
  for (b_k = 0; b_k < 49; b_k++) {
    d[b_k] += S_in[b_k];
  }
  st.site = &e_emlrtRSI;
  c_y = det(&st, d);
  st.site = &e_emlrtRSI;
  if (c_y < 0.0) {
    emlrtErrorWithMessageIdR2018a(
        &st, &b_emlrtRTEI, "Coder:toolbox:ElFunDomainError",
        "Coder:toolbox:ElFunDomainError", 3, 4, 4, "sqrt");
  }
  c_y = muDoubleScalarSqrt(c_y);
  a = sigma2f / c_y * var_short;
  for (b_k = 0; b_k < 140; b_k++) {
    T[b_k] *= inp[b_k];
  }
  sum(T, x);
  for (k = 0; k < 20; k++) {
    x[k] = muDoubleScalarExp(-x[k] / 2.0);
  }
  /*  short-term mean */
  for (k = 0; k < 7; k++) {
    for (b_k = 0; b_k < 20; b_k++) {
      inp_tmp = b_k + 20 * k;
      inp_s[inp_tmp] = Smodel->Xu[inp_tmp] - m_in[k];
    }
  }
  for (b_k = 0; b_k < 20; b_k++) {
    c_y = 0.0;
    for (i = 0; i < 20; i++) {
      c_y += Smodel->invKm[b_k + 20 * i] * Smodel->postm[i];
    }
    beta_s[b_k] = c_y;
  }
  for (k = 0; k < 7; k++) {
    d_y[k] = muDoubleScalarExp(-2.0 * Smodel->GP.logtheta[k]);
  }
  memset(&iLambda_s[0], 0, 49U * sizeof(real_T));
  /*  DxD, inversed squared length scale */
  for (inp_tmp = 0; inp_tmp < 7; inp_tmp++) {
    iLambda_s[inp_tmp + 7 * inp_tmp] = d_y[inp_tmp];
    d_y[inp_tmp] = muDoubleScalarExp(2.0 * Smodel->GP.logtheta[inp_tmp]);
  }
  memset(&d[0], 0, 49U * sizeof(real_T));
  for (inp_tmp = 0; inp_tmp < 7; inp_tmp++) {
    d[inp_tmp + 7 * inp_tmp] = d_y[inp_tmp];
    for (b_k = 0; b_k < 7; b_k++) {
      c_y = 0.0;
      for (i = 0; i < 7; i++) {
        c_y += S_in[inp_tmp + 7 * i] * iLambda_s[i + 7 * b_k];
      }
      y_tmp[inp_tmp + 7 * b_k] = c_y;
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
  st.site = &f_emlrtRSI;
  mldivide(&st, c_I, y_tmp);
  for (b_k = 0; b_k < 49; b_k++) {
    y_tmp[b_k] = (real_T)iR_tmp[b_k] - y_tmp[b_k];
  }
  st.site = &g_emlrtRSI;
  for (b_k = 0; b_k < 7; b_k++) {
    for (i = 0; i < 7; i++) {
      c_y = 0.0;
      for (inp_tmp = 0; inp_tmp < 7; inp_tmp++) {
        c_y += iLambda_s[b_k + 7 * inp_tmp] * y_tmp[inp_tmp + 7 * i];
      }
      c_I[b_k + 7 * i] = c_y;
    }
  }
  b_mtimes(inp_s, c_I, T);
  st.site = &h_emlrtRSI;
  var_short = b_sumColumnB(*(real_T(*)[7]) & Smodel->GP.logtheta[0]);
  var_short = muDoubleScalarExp(var_short);
  for (b_k = 0; b_k < 49; b_k++) {
    d[b_k] += S_in[b_k];
  }
  st.site = &h_emlrtRSI;
  c_y = det(&st, d);
  st.site = &h_emlrtRSI;
  if (c_y < 0.0) {
    emlrtErrorWithMessageIdR2018a(
        &st, &b_emlrtRTEI, "Coder:toolbox:ElFunDomainError",
        "Coder:toolbox:ElFunDomainError", 3, 4, 4, "sqrt");
  }
  c_y = muDoubleScalarSqrt(c_y);
  var_short *= sigma2f_s / c_y;
  for (b_k = 0; b_k < 140; b_k++) {
    T[b_k] *= inp_s[b_k];
  }
  sum(T, b_x);
  for (k = 0; k < 20; k++) {
    b_x[k] = muDoubleScalarExp(-b_x[k] / 2.0);
  }
  /*  DGP predict mean */
  st.site = &i_emlrtRSI;
  for (b_k = 0; b_k < 20; b_k++) {
    x[b_k] = a * x[b_k] * beta[b_k];
  }
  c_y = c_sumColumnB(x);
  st.site = &i_emlrtRSI;
  for (b_k = 0; b_k < 20; b_k++) {
    b_x[b_k] = var_short * b_x[b_k] * beta_s[b_k];
  }
  var_short = c_sumColumnB(b_x);
  *mu_MM = c_y + var_short;
  /*  long-term  var */
  for (k = 0; k < 7; k++) {
    d_y[k] = muDoubleScalarExp(model->GP.logtheta[k]);
    for (b_k = 0; b_k < 20; b_k++) {
      inp_tmp = b_k + 20 * k;
      T[inp_tmp] = inp[inp_tmp] / d_y[k];
    }
  }
  /*  xu-mu/sqrt(iLambda) */
  for (b_k = 0; b_k < 140; b_k++) {
    c_y = T[b_k];
    c_y *= c_y;
    T[b_k] = c_y;
  }
  sum(T, log_k);
  c_y = 2.0 * model->GP.logtheta[7];
  for (b_k = 0; b_k < 20; b_k++) {
    log_k[b_k] = c_y - log_k[b_k] / 2.0;
  }
  /*  DxD */
  for (k = 0; k < 7; k++) {
    c_y = model->GP.logtheta[k];
    d_y[k] = muDoubleScalarExp(2.0 * c_y);
    for (b_k = 0; b_k < 20; b_k++) {
      inp_tmp = b_k + 20 * k;
      Zeta_i[inp_tmp] = inp[inp_tmp] / d_y[k];
    }
    var_short = muDoubleScalarExp(-2.0 * c_y);
    d_y[k] = var_short + var_short;
  }
  memset(&y_tmp[0], 0, 49U * sizeof(real_T));
  for (inp_tmp = 0; inp_tmp < 7; inp_tmp++) {
    y_tmp[inp_tmp + 7 * inp_tmp] = d_y[inp_tmp];
  }
  for (b_k = 0; b_k < 7; b_k++) {
    for (i = 0; i < 7; i++) {
      c_y = 0.0;
      for (inp_tmp = 0; inp_tmp < 7; inp_tmp++) {
        c_y += S_in[b_k + 7 * inp_tmp] * y_tmp[inp_tmp + 7 * i];
      }
      inp_tmp = b_k + 7 * i;
      iLambda[inp_tmp] = c_y + (real_T)iR_tmp[inp_tmp];
    }
  }
  st.site = &j_emlrtRSI;
  c_y = det(&st, iLambda);
  st.site = &j_emlrtRSI;
  if (c_y < 0.0) {
    emlrtErrorWithMessageIdR2018a(
        &st, &b_emlrtRTEI, "Coder:toolbox:ElFunDomainError",
        "Coder:toolbox:ElFunDomainError", 3, 4, 4, "sqrt");
  }
  c_y = muDoubleScalarSqrt(c_y);
  a = 1.0 / c_y;
  st.site = &k_emlrtRSI;
  for (b_k = 0; b_k < 140; b_k++) {
    T[b_k] = -Zeta_i[b_k];
  }
  memcpy(&iLambda_s[0], &S_in[0], 49U * sizeof(real_T));
  b_st.site = &k_emlrtRSI;
  mldivide(&b_st, iLambda, iLambda_s);
  for (b_k = 0; b_k < 49; b_k++) {
    iLambda_s[b_k] /= 2.0;
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
  b_mtimes(Zeta_i, iLambda_s, aQ);
  b_mtimes(T, iLambda_s, e_y);
  for (b_k = 0; b_k < 140; b_k++) {
    dv[b_k] = 2.0 * aQ[b_k];
  }
  c_mtimes(dv, T, A);
  for (k = 0; k < 20; k++) {
    for (b_k = 0; b_k < 20; b_k++) {
      c[b_k + 20 * k] = log_k[b_k] + log_k[k];
    }
  }
  for (b_k = 0; b_k < 140; b_k++) {
    aQ[b_k] *= Zeta_i[b_k];
    e_y[b_k] *= T[b_k];
  }
  sum(aQ, x);
  sum(e_y, b_x);
  for (k = 0; k < 20; k++) {
    for (b_k = 0; b_k < 20; b_k++) {
      b_c[b_k + 20 * k] = x[b_k] + b_x[k];
    }
  }
  for (k = 0; k < 400; k++) {
    c[k] = muDoubleScalarExp(c[k] + (b_c[k] - A[k]));
  }
  for (b_k = 0; b_k < 20; b_k++) {
    for (i = 0; i < 20; i++) {
      inp_tmp = i + 20 * b_k;
      A[inp_tmp] = beta[i] * beta[b_k] - (model->invKm[inp_tmp] - y[inp_tmp]);
    }
  }
  for (b_k = 0; b_k < 400; b_k++) {
    A[b_k] *= a * c[b_k];
  }
  st.site = &l_emlrtRSI;
  for (inp_tmp = 0; inp_tmp < 20; inp_tmp++) {
    x[inp_tmp] = sumColumnB(A, inp_tmp + 1);
  }
  st.site = &l_emlrtRSI;
  var_long = c_sumColumnB(x);
  var_long += sigma2f;
  /*  short-term var */
  for (k = 0; k < 7; k++) {
    d_y[k] = muDoubleScalarExp(Smodel->GP.logtheta[k]);
    for (b_k = 0; b_k < 20; b_k++) {
      inp_tmp = b_k + 20 * k;
      T[inp_tmp] = inp_s[inp_tmp] / d_y[k];
    }
  }
  /*  xu-mu/sqrt(iLambda) */
  for (b_k = 0; b_k < 140; b_k++) {
    c_y = T[b_k];
    c_y *= c_y;
    T[b_k] = c_y;
  }
  sum(T, log_ks);
  c_y = 2.0 * Smodel->GP.logtheta[7];
  for (b_k = 0; b_k < 20; b_k++) {
    log_ks[b_k] = c_y - log_ks[b_k] / 2.0;
  }
  /*  DxD */
  for (k = 0; k < 7; k++) {
    c_y = Smodel->GP.logtheta[k];
    d_y[k] = muDoubleScalarExp(2.0 * c_y);
    for (b_k = 0; b_k < 20; b_k++) {
      inp_tmp = b_k + 20 * k;
      inp[inp_tmp] = inp_s[inp_tmp] / d_y[k];
    }
    var_short = muDoubleScalarExp(-2.0 * c_y);
    d_y[k] = var_short + var_short;
  }
  memset(&y_tmp[0], 0, 49U * sizeof(real_T));
  for (inp_tmp = 0; inp_tmp < 7; inp_tmp++) {
    y_tmp[inp_tmp + 7 * inp_tmp] = d_y[inp_tmp];
  }
  for (b_k = 0; b_k < 7; b_k++) {
    for (i = 0; i < 7; i++) {
      c_y = 0.0;
      for (inp_tmp = 0; inp_tmp < 7; inp_tmp++) {
        c_y += S_in[b_k + 7 * inp_tmp] * y_tmp[inp_tmp + 7 * i];
      }
      inp_tmp = b_k + 7 * i;
      iLambda[inp_tmp] = c_y + (real_T)iR_tmp[inp_tmp];
    }
  }
  st.site = &m_emlrtRSI;
  c_y = det(&st, iLambda);
  st.site = &m_emlrtRSI;
  if (c_y < 0.0) {
    emlrtErrorWithMessageIdR2018a(
        &st, &b_emlrtRTEI, "Coder:toolbox:ElFunDomainError",
        "Coder:toolbox:ElFunDomainError", 3, 4, 4, "sqrt");
  }
  c_y = muDoubleScalarSqrt(c_y);
  a = 1.0 / c_y;
  st.site = &n_emlrtRSI;
  for (b_k = 0; b_k < 140; b_k++) {
    T[b_k] = -inp[b_k];
  }
  memcpy(&iLambda_s[0], &S_in[0], 49U * sizeof(real_T));
  b_st.site = &n_emlrtRSI;
  mldivide(&b_st, iLambda, iLambda_s);
  for (b_k = 0; b_k < 49; b_k++) {
    iLambda_s[b_k] /= 2.0;
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
  b_mtimes(inp, iLambda_s, aQ);
  b_mtimes(T, iLambda_s, e_y);
  for (b_k = 0; b_k < 140; b_k++) {
    dv[b_k] = 2.0 * aQ[b_k];
  }
  c_mtimes(dv, T, A);
  for (k = 0; k < 20; k++) {
    for (b_k = 0; b_k < 20; b_k++) {
      c[b_k + 20 * k] = log_ks[b_k] + log_ks[k];
    }
  }
  for (b_k = 0; b_k < 140; b_k++) {
    aQ[b_k] *= inp[b_k];
    e_y[b_k] *= T[b_k];
  }
  sum(aQ, x);
  sum(e_y, b_x);
  for (k = 0; k < 20; k++) {
    for (b_k = 0; b_k < 20; b_k++) {
      b_c[b_k + 20 * k] = x[b_k] + b_x[k];
    }
  }
  for (k = 0; k < 400; k++) {
    c[k] = muDoubleScalarExp(c[k] + (b_c[k] - A[k]));
  }
  for (b_k = 0; b_k < 20; b_k++) {
    for (i = 0; i < 20; i++) {
      inp_tmp = i + 20 * b_k;
      A[inp_tmp] =
          beta_s[i] * beta_s[b_k] - (Smodel->invKm[inp_tmp] - b_y[inp_tmp]);
    }
  }
  for (b_k = 0; b_k < 400; b_k++) {
    A[b_k] *= a * c[b_k];
  }
  st.site = &o_emlrtRSI;
  for (inp_tmp = 0; inp_tmp < 20; inp_tmp++) {
    x[inp_tmp] = sumColumnB(A, inp_tmp + 1);
  }
  st.site = &o_emlrtRSI;
  var_short = c_sumColumnB(x);
  if (var_short > sigma2f_s) {
    var_short = sigma2f_s;
  }
  var_short += sigma2f_s;
  /*  cross var */
  for (k = 0; k < 7; k++) {
    d_y[k] = muDoubleScalarExp(-2.0 * model->GP.logtheta[k]) +
             muDoubleScalarExp(-2.0 * Smodel->GP.logtheta[k]);
  }
  memset(&y_tmp[0], 0, 49U * sizeof(real_T));
  for (inp_tmp = 0; inp_tmp < 7; inp_tmp++) {
    y_tmp[inp_tmp + 7 * inp_tmp] = d_y[inp_tmp];
  }
  for (b_k = 0; b_k < 7; b_k++) {
    for (i = 0; i < 7; i++) {
      c_y = 0.0;
      for (inp_tmp = 0; inp_tmp < 7; inp_tmp++) {
        c_y += S_in[b_k + 7 * inp_tmp] * y_tmp[inp_tmp + 7 * i];
      }
      inp_tmp = b_k + 7 * i;
      iLambda[inp_tmp] = c_y + (real_T)iR_tmp[inp_tmp];
    }
  }
  st.site = &p_emlrtRSI;
  c_y = det(&st, iLambda);
  st.site = &p_emlrtRSI;
  if (c_y < 0.0) {
    emlrtErrorWithMessageIdR2018a(
        &st, &b_emlrtRTEI, "Coder:toolbox:ElFunDomainError",
        "Coder:toolbox:ElFunDomainError", 3, 4, 4, "sqrt");
  }
  c_y = muDoubleScalarSqrt(c_y);
  a = 1.0 / c_y;
  st.site = &q_emlrtRSI;
  memcpy(&iLambda_s[0], &S_in[0], 49U * sizeof(real_T));
  b_st.site = &q_emlrtRSI;
  mldivide(&b_st, iLambda, iLambda_s);
  for (b_k = 0; b_k < 49; b_k++) {
    iLambda_s[b_k] /= 2.0;
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
  b_mtimes(Zeta_i, iLambda_s, aQ);
  b_mtimes(T, iLambda_s, e_y);
  for (b_k = 0; b_k < 140; b_k++) {
    dv[b_k] = 2.0 * aQ[b_k];
  }
  c_mtimes(dv, T, A);
  for (k = 0; k < 20; k++) {
    for (b_k = 0; b_k < 20; b_k++) {
      c[b_k + 20 * k] = log_k[b_k] + log_ks[k];
    }
  }
  for (b_k = 0; b_k < 140; b_k++) {
    aQ[b_k] *= Zeta_i[b_k];
    e_y[b_k] *= T[b_k];
  }
  sum(aQ, x);
  sum(e_y, b_x);
  for (k = 0; k < 20; k++) {
    for (b_k = 0; b_k < 20; b_k++) {
      b_c[b_k + 20 * k] = x[b_k] + b_x[k];
    }
  }
  for (k = 0; k < 400; k++) {
    c[k] = muDoubleScalarExp(c[k] + (b_c[k] - A[k]));
  }
  for (b_k = 0; b_k < 20; b_k++) {
    for (i = 0; i < 20; i++) {
      A[i + 20 * b_k] = beta[i] * beta_s[b_k];
    }
  }
  for (b_k = 0; b_k < 400; b_k++) {
    A[b_k] *= a * c[b_k];
  }
  st.site = &r_emlrtRSI;
  for (inp_tmp = 0; inp_tmp < 20; inp_tmp++) {
    x[inp_tmp] = sumColumnB(A, inp_tmp + 1);
  }
  st.site = &r_emlrtRSI;
  c_y = c_sumColumnB(x);
  /*  total var = long + short + cross */
  *var_MM = ((var_long + var_short) + 2.0 * c_y) - *mu_MM * *mu_MM;
  if (*var_MM < 0.0) {
    *var_MM = 0.0;
  }
  /*  in case for som numerical instabilities */
}

/* End of code generation (DGP_varsgpPredict_y_c.c) */
