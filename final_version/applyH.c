//Update A = H*A, where H is a Householder matrix with a Householder vector v normalized
//to have v(1) = 1. 

#include <stdio.h>
#include <stdlib.h>
#include "cblas.h"
#include <math.h>

double *applyH(double *A, int mA, int nA, double *v, int mv, int nv)
{

        double *w, *u;
        double gamma;
        double alpha;
        double beta;
        int i, row, col;
        int inc = 1;

	if (mv != mA || nv > 1)
	{
		printf("Size is wrong, please reset the input matrix A and vector v.\n");
		double *p = NULL;
		A = p;
	}

	//Set up a 1*nA vector w and initialize it with all 0s.
	w = (double *) malloc(1*nA*sizeof(double));
	for(i=0;i<nA;i++)
	{
		w[i] = 0;
	}

	//Set up a mv*1 vector u and initialize it with all 0s.
	u = (double *) malloc(mv*1*sizeof(double));
	for(i=0;i<mv;i++)
	{
		u[i] = 0;
	}

	//HOW TO IMPLEMENT THIS WITH BLAS FUNCTION
	//Original computation was one line:
	//
	//      A = A - (2/(v'*v))*v*(v'*A);
	//
	//The intermediate temporary 1D array w is needed in C++/Fortran,
	//but the vector v is not. You can use a temporary arrav v if it
	//makes things easier (which it will).
	//
	//gamma = 2/(v'*v);  % dotproduct of v with itself (don't use dnrm2() squared)
	//w     = v'*A;          % call to dgemv()
	//u     = gamma*v;       % call to dscal()
	//A     = A - u*w;       % call to dger()

        gamma = 2/(cblas_ddot(mv, v, inc, v, inc));
        alpha= 1.0;
        beta= 0.0;

	//Transpose the matrix.
	//When I used the original matrix and set the argument to CblasTrans, the result is wrong.
	//So here, I transpose in advance and send this transposed matrix to cblas_dgemv.
	double *A_ = (double *) malloc(nA*mA*sizeof(double));
	for(row=0;row<nA;row++)
	{
		for(col=0;col<mA;col++)
		{
			A_[row*mA+col] = A[col*nA+row];
		}
	}

	//Multiply transposed matrix A with vector v and return the result to w.
        cblas_dgemv(CblasRowMajor, CblasNoTrans, nA, mA, alpha, A_, mA, v, inc, beta, w, inc);

	//Multiply gamma and vector v, store the result to vector u.
        cblas_dscal(mv, gamma, v, inc);
        
        for(i=0;i<mv;i++)
        {
            u[i] = v[i];
        }

	//Use cblas_dger to compute A = alpha*u*w'+A, i.e., A = A-u*w.
        alpha = -1.0;
        cblas_dger(CblasRowMajor, mA, nA, alpha, u, inc, w, inc, A, nA);

        return A;
}
