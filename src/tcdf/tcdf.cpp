/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: main.c
 *
 * MATLAB Coder version            : 3.3
 * C/C++ source code generated on  : 24-Apr-2018 18:23:19
 */

/*************************************************************************/
/* This automatically generated example C main file shows how to call    */
/* entry-point functions that MATLAB Coder generated. You must customize */
/* this file for your application. Do not modify this file directly.     */
/* Instead, make a copy of this file, modify it, and integrate it into   */
/* your development environment.                                         */
/*                                                                       */
/* This file initializes entry-point function arguments to a default     */
/* size and value before calling the entry-point functions. It does      */
/* not store or use any values returned from the entry-point functions.  */
/* If necessary, it does pre-allocate memory for returned values.        */
/* You can use this file as a starting point for a main function that    */
/* you can deploy in your application.                                   */
/*                                                                       */
/* After you copy the file, and before you deploy it, you must make the  */
/* following changes:                                                    */
/* * For variable-size function arguments, change the example sizes to   */
/* the sizes that your application requires.                             */
/* * Change the example values of function arguments to the values that  */
/* your application requires.                                            */
/* * If the entry-point functions return values, store these values or   */
/* otherwise use them as required by your application.                   */
/*                                                                       */
/*************************************************************************/
/* Include Files */
#include "tcdf.h"
#include "ctcdf.h"

#include "rt_nonfinite.h"
#include "ctcdf_terminate.h"
#include "ctcdf_emxAPI.h"
#include "ctcdf_initialize.h"

namespace matlabautogen
{


/*
 * Arguments    : void
 * Return Type  : emxArray_real_T *
 */
emxArray_real_T *c_argInit_UnboundedxUnbounded_r(const Eigen::MatrixXd &mat)
{
  emxArray_real_T *result;
  int iv0[2] = { (int)mat.rows(), (int)mat.cols() };

 

  /* Set the size of the array.
     Change this size to the value that the application requires. */
  result = emxCreateND_real_T(2, iv0);
  memcpy(result->data, mat.data(), mat.size());

 // int idx0;
 //  int idx1;
 //    /* Loop over the array to initialize each element. */
 //  for (idx0 = 0; idx0 < result->size[0U]; idx0++) {
 //    for (idx1 = 0; idx1 < result->size[1U]; idx1++) {
 //      /* Set the value of the array element.
 //         Change this value to the value that the application requires. */
 //      result->data[idx0 + result->size[0] * idx1] = argInit_real_T();
 //    }
 //  }

  return result;
}

/*
 * Arguments    : void
 * Return Type  : void
 */
void main_ctcdf(const Eigen::MatrixXd &x_eigen, const Eigen::MatrixXd &v_eigen, Eigen::MatrixXd &res_eigen)
{
  emxArray_real_T *res;
  emxArray_real_T *x;
  emxArray_real_T *v;
  emxInitArray_real_T(&res, 2);

  /* Initialize function 'ctcdf' input arguments. */
  /* Initialize function input argument 'x'. */
  x = c_argInit_UnboundedxUnbounded_r(x_eigen);

  /* Initialize function input argument 'v'. */
  v = c_argInit_UnboundedxUnbounded_r(v_eigen);

  /* Call the entry-point 'ctcdf'. */
  ctcdf(x, v, res);

  res_eigen = Eigen::Map<Eigen::MatrixXd>(res->data, res->size[0], res->size[1]);

  emxDestroyArray_real_T(res);
  emxDestroyArray_real_T(v);
  emxDestroyArray_real_T(x);
}

/*
 * Arguments    : int argc
 *                const char * const argv[]
 * Return Type  : int
 */
void tcdf(const Eigen::MatrixXd &x, const Eigen::MatrixXd &v, Eigen::MatrixXd &res)
{
  /* Initialize the application.
     You do not need to do this more than one time. */
  ctcdf_initialize();

  /* Invoke the entry-point functions.
     You can call entry-point functions multiple times. */
  main_ctcdf(x, v, res);

  /* Terminate the application.
     You do not need to do this more than one time. */
  ctcdf_terminate();
}
}
/*
 * File trailer for main.c
 *
 * [EOF]
 */
