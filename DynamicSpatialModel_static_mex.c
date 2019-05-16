/*
 * DynamicSpatialModel_static_mex.c : The gateway routine used to call the Static function located in DynamicSpatialModel_static.c
 *
 * Warning : this file is generated automatically by Dynare
 *           from model file (.mod)

 */

#include "mex.h"

void Static(double *y, double *x, int nb_row_x, double *params, double *residual, double *g1, double *v2);
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  double *y, *x, *params;
  double *residual, *g1, *v2;
  int nb_row_x;

  /* Create a pointer to the input matrix y. */
  y = mxGetPr(prhs[0]);

  /* Create a pointer to the input matrix x. */
  x = mxGetPr(prhs[1]);

  /* Create a pointer to the input matrix params. */
  params = mxGetPr(prhs[2]);

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
      plhs[1] = mxCreateDoubleMatrix(262, 262, mxREAL);
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

  /* Call the C subroutines. */
  Static(y, x, nb_row_x, params, residual, g1, v2);
}

