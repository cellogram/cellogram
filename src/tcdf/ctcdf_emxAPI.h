/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: ctcdf_emxAPI.h
 *
 * MATLAB Coder version            : 3.3
 * C/C++ source code generated on  : 24-Apr-2018 18:23:19
 */

#ifndef CTCDF_EMXAPI_H
#define CTCDF_EMXAPI_H

/* Include Files */
#include <math.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include "rt_nonfinite.h"
#include "rtwtypes.h"
#include "ctcdf_types.h"
namespace matlabautogen
{
/* Function Declarations */
 emxArray_real_T *emxCreateND_real_T(int numDimensions, int *size);
 emxArray_real_T *emxCreateWrapperND_real_T(double *data, int
  numDimensions, int *size);
 emxArray_real_T *emxCreateWrapper_real_T(double *data, int rows, int cols);
 emxArray_real_T *emxCreate_real_T(int rows, int cols);
 void emxDestroyArray_real_T(emxArray_real_T *emxArray);
 void emxInitArray_real_T(emxArray_real_T **pEmxArray, int numDimensions);
}
#endif

/*
 * File trailer for ctcdf_emxAPI.h
 *
 * [EOF]
 */
