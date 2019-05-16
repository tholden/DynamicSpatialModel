/*
 * DynamicSpatialModel_dynamic_mex.c : The gateway routine used to call the Dynamic function located in DynamicSpatialModel_dynamic.c
 *
 * Warning : this file is generated automatically by Dynare
 *           from model file (.mod)

 */

#include "mex.h"

void Dynamic(double *y, double *x, int nb_row_x, double *params, double *steady_state, int it_, double *residual, double *g1, double *v2, double *v3);
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  double *y, *x, *params, *steady_state;
  double *residual, *g1, *v2, *v3;
  int nb_row_x, it_;

  /* Check that no derivatives of higher order than computed are being requested */
  if (nlhs > 2)
    mexErrMsgTxt("Derivatives of higher order than computed have been requested");
  /* Create a pointer to the input matrix y. */
  y = mxGetPr(prhs[0]);

  /* Create a pointer to the input matrix x. */
  x = mxGetPr(prhs[1]);

  /* Create a pointer to the input matrix params. */
  params = mxGetPr(prhs[2]);

  /* Create a pointer to the input matrix steady_state. */
  steady_state = mxGetPr(prhs[3]);

  /* Fetch time index */
  it_ = (int) mxGetScalar(prhs[4]) - 1;

  /* Gets number of rows of matrix x. */
  nb_row_x = mxGetM(prhs[1]);

  residual = NULL;
  if (nlhs >= 1)
  {
     /* Set the output pointer to the output matrix residual. */
     plhs[0] = mxCreateDoubleMatrix(262,1, mxREAL);
     /* Create a C pointer to a copy of the output matrix residual. */
     residual = mxGetPr(plhs[0]);
  }

  g1 = NULL;
  if (nlhs >= 2)
  {
     /* Set the output pointer to the output matrix g1. */
     plhs[1] = mxCreateDoubleMatrix(262, 538, mxREAL);
     /* Create a C pointer to a copy of the output matrix g1. */
     g1 = mxGetPr(plhs[1]);
  }

  v2 = NULL;
 if (nlhs >= 3)
  {
     /* Set the output pointer to the output matrix v2. */
     plhs[2] = mxCreateDoubleMatrix(0, 3, mxREAL);
     v2 = mxGetPr(plhs[2]);
  }

  v3 = NULL;
 if (nlhs >= 4)
  {
     /* Set the output pointer to the output matrix v3. */
     plhs[3] = mxCreateDoubleMatrix(0, 3, mxREAL);
     v3 = mxGetPr(plhs[3]);
  }

  /* Call the C subroutines. */
  Dynamic(y, x, nb_row_x, params, steady_state, it_, residual, g1, v2, v3);
}
