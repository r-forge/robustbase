#include "mex.h"
#include "mlmc.c"
#include<stdio.h>

void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray*prhs[] )

{
    double *yout;
    double *yin;
    long m,n,i;

    /* Check for proper number of arguments */

    if (nrhs != 1) {
        mexErrMsgTxt("One input argument required.");
    } else if (nlhs > 1) {
        mexErrMsgTxt("Too many output arguments.");
    }

    m = mxGetM(prhs[0]);
    n = mxGetN(prhs[0]);

    if (n!=1 && m==1)
    {
        mexErrMsgTxt("Input must be a columnvector.");
    }
    /* Create a matrix for the return argument */
    plhs[0] = mxCreateDoubleMatrix(1, n, mxREAL);

    /* Assign pointers to the various parameters */
    yout = mxGetPr(plhs[0]);
    yin = mxGetPr(prhs[0]);

        for (i=0;i<n;i++)
                mlmc(&yout[i],&yin[i*m],&m);

    return;

}
