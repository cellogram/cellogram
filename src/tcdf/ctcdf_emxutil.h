/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: ctcdf_emxutil.h
 *
 * MATLAB Coder version            : 3.3
 * C/C++ source code generated on  : 24-Apr-2018 18:23:19
 */

#ifndef CTCDF_EMXUTIL_H
#define CTCDF_EMXUTIL_H

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
 void emxEnsureCapacity(matlabautogen::emxArray__common *emxArray, int oldNumel, unsigned
  int elementSize);
 void emxFree_real_T(matlabautogen::emxArray_real_T **pEmxArray);
 void emxInit_real_T(matlabautogen::emxArray_real_T **pEmxArray, int numDimensions);
}
#endif

/*
 * File trailer for ctcdf_emxutil.h
 *
 * [EOF]
 */
