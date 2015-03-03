//Update A = H*A, where H is a Householder matrix with a Householder vector v normalized
//to have v(1) = 1. 

#include <stdio.h>
#include <stdlib.h>
#include "cblas.h"
#include <math.h>

double *applyH(double *A, int mA, int nA, double *v, int mv, int nv)
{

	printf("~~~~~~~~~~applyH function start~~~~~~~~~~\n");

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
	printf("Print gamma = %lf\n", gamma);
        alpha= 1.0;
        beta= 0.0;

	//Print out the input matrix A and input vector v.
	printf("Let's test the input matrix A and input vector v.\n");
	for(row=0;row<mA;row++)
	{
		for(col=0;col<nA;col++)
		{
			printf("A[%d][%d] = %lf ", row, col, A[row*nA+col]);
		}
		printf("\n");
	}
	for(i=0;i<mv;i++)
		printf("v[%d] = %lf ", i, v[i]);
	printf("\n");
	printf("alpha = %lf, beta = %lf, inc = %d\n", alpha, beta, inc);

	//Transpose the matrix and print it out.
	//When I used the original matrix and set the argument to CblasTrans, the result is wrong.
	//So here, I transpose in advance and send this transposed matrix to cblas_dgemv.
	double *A_ = (double *) malloc(nA*mA*sizeof(double));
	printf("Let's test the transposed A, A'.\n");
	for(row=0;row<nA;row++)
	{
		for(col=0;col<mA;col++)
		{
			A_[row*mA+col] = A[col*nA+row];
			printf("A'[%d][%d] = %lf ", row, col, A_[row*mA+col]);
		}
		printf("\n");
	}

	//Multiply transposed matrix A with vector v and return the result to w.
        cblas_dgemv(CblasRowMajor, CblasNoTrans, nA, mA, alpha, A_, mA, v, inc, beta, w, inc);

	/*
	//Test the cblas_dgemv function alone
	//double A2[12] = {1.0, 1.0, 1.0, 1.0, -1.0, 4.0, 4.0, -1.0, 4.0, -2.0, 2.0, 0.0};
	double A2[12] = {1.0, -1.0, 4.0, 1.0, 4.0, -2.0, 1.0, 4.0, 2.0, 1.0, -1.0, 0.0};
	double v2[4] = {-1.0, 1.0, 1.0, 1.0};
	double w2[3] = {0.0, 0.0, 0.0};
	cblas_dgemv(CblasRowMajor, CblasTrans, 3, 4, alpha, A2, 1, v2, inc, beta, w2, inc);
	printf("Test the cblas_dgemv function alone.\n");
	int n2;
	for(n2=0;n2<3;n2++)
		printf("w2[%d] = %lf\t", n2, w2[n2]);
	printf("\n");
	*/
	
	//Print out the cblas_dgemv function result w.
	printf("Let's test w.\n");
	for(i=0;i<nA;i++)
		printf("w[%d] = %lf ", i, w[i]);
	printf("\n");

	//Multiply gamma and vector v, store the result to vector u, and print it out.
        cblas_dscal(mv, gamma, v, inc);
        printf("Let's test u.\n");
        for(i=0;i<mv;i++)
        {
            u[i] = v[i];
            printf("u[%d] = %lf ", i, u[i]);
        }
        printf("\n");

	printf("Before update A with \"A = A - u * w\", let's print the original A.\n");
	for(row=0;row<mA;row++)
	{
		for(col=0;col<nA;col++)
		{
			printf("A[%d][%d] = %lf ", row, col, A[row*nA+col]);
		}
		printf("\n");
	}

	//Use cblas_dger to compute A = alpha*u*w'+A, i.e., A = A-u*w.
        alpha = -1.0;
        cblas_dger(CblasRowMajor, mA, nA, alpha, u, inc, w, inc, A, nA);
        
        /*
        printf("Test the cblas_dger function alone.\n");
        double test_A[6] = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0};
        double test_u[3] = {1.0, 1.0, 2.0};
        double test_w[2] = {2.0, 3.0};
        cblas_dger(CblasRowMajor, 3, 2, -1.0, test_u, 1, test_w, 1, test_A, 2);
        for(row=0;row<3;row++)
	{
		for(col=0;col<2;col++)
		{
			printf("test_A[%d][%d] = %lf\t", row, col, test_A[row*2+col]);
		}
		printf("\n");
	}
	*/

	printf("Print the result of A.\n");
	for(row=0;row<mA;row++)
	{
		for(col=0;col<nA;col++)
		{
			printf("A[%d][%d] = %lf ", row, col, A[row*nA+col]);
		}
		printf("\n");
	}
	
	printf("~~~~~~~~~~applyH function end~~~~~~~~~~\n");

        return A;
}
