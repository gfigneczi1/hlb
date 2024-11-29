/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * DGP_varsgpPredict_y_c_types.h
 *
 * Code generation for function 'DGP_varsgpPredict_y_c'
 *
 */

#pragma once

/* Include files */
#include "rtwtypes.h"
#include "emlrt.h"

/* Type Definitions */
#ifndef typedef_struct1_T
#define typedef_struct1_T
typedef struct {
  char_T type[8];
  real_T nParams;
  real_T logtheta;
} struct1_T;
#endif /* typedef_struct1_T */

#ifndef typedef_struct2_T
#define typedef_struct2_T
typedef struct {
  char_T type[5];
  real_T logtheta[8];
  real_T nParams;
  real_T constDiag;
} struct2_T;
#endif /* typedef_struct2_T */

#ifndef typedef_struct0_T
#define typedef_struct0_T
typedef struct {
  char_T type[7];
  char_T objectFunc[3];
  struct1_T Likelihood;
  struct2_T GP;
  real_T n;
  real_T D;
  real_T y[996];
  real_T X[6972];
  real_T m;
  char_T indType[9];
  char_T indRepar[2];
  real_T Xu[140];
  real_T Xuinit[140];
  real_T nIndParams;
  real_T yy;
  real_T vary;
  real_T jitter;
  real_T sigma2;
  real_T postm[20];
  real_T postS[400];
  real_T alpha[20];
  real_T Kmm[400];
  real_T Kmn[19920];
  real_T Knm[19920];
  real_T L[400];
  real_T invKm[400];
} struct0_T;
#endif /* typedef_struct0_T */

#ifndef typedef_struct3_T
#define typedef_struct3_T
typedef struct {
  char_T type[7];
  char_T objectFunc[3];
  struct1_T Likelihood;
  struct2_T GP;
  real_T n;
  real_T D;
  real_T y[996];
  real_T X[6972];
  real_T m;
  char_T indType[9];
  char_T indRepar[2];
  real_T Xu[140];
  real_T Xuinit[140];
  real_T nIndParams;
  real_T yy;
  real_T vary;
  real_T jitter;
  real_T postm[20];
  real_T postS[400];
  real_T Kmm[400];
  real_T Kmn[19920];
  real_T Knm[19920];
  real_T L[400];
  real_T invKm[400];
} struct3_T;
#endif /* typedef_struct3_T */

#ifndef typedef_b_DGP_varsgpPredict_y_c_api
#define typedef_b_DGP_varsgpPredict_y_c_api
typedef struct {
  struct0_T model;
  struct3_T Smodel;
} b_DGP_varsgpPredict_y_c_api;
#endif /* typedef_b_DGP_varsgpPredict_y_c_api */

#ifndef typedef_DGP_varsgpPredict_y_cStackData
#define typedef_DGP_varsgpPredict_y_cStackData
typedef struct {
  b_DGP_varsgpPredict_y_c_api f0;
} DGP_varsgpPredict_y_cStackData;
#endif /* typedef_DGP_varsgpPredict_y_cStackData */

/* End of code generation (DGP_varsgpPredict_y_c_types.h) */
