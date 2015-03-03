//Perform a complete QR factorization on an mxn matrix A, returning the upper
//triangular matrix R in the upper triangular part of the array A, and storing the
//Householder vectors v_k in the strictly lower triangular part of the array A.
//The first entry of each Householder vector is not stored since it is assumed
//that v_k(1) = 1 for all of the Householder vectors v_k.

#include <stdio.h>
#include <stdlib.h>
#include "cblas.h"
#include "house.h"
#include "applyH.h"
#include <math.h>

double *qrhouse(double *A, int m, int n)
{
	int i;
	int row;
	int col;
	int inc = 1;
	double alpha;
	double beta;
	
	//Create an m size vector v and initialize it with all 0s.
	double *v = (double *) malloc(m*1*sizeof(double));
	for(i=0;i<m;i++)
		v[i] = 0;
	
	int min;
	if(m-1<n)
		min = m-1;
	else
		min = n;
	int k;
	for(k=0;k<min;k++)
	{
		//Get the first row of subset of A.
		double *sub_Av = (double *) malloc((m-k)*1*sizeof(double));
		for(i=0;i<m-k;i++)
		{
			sub_Av[i] = A[(i+k)*n+k];
		}
		
		//If only the first element of this vector is not 0, the reflector will be all 0s and thus H will be I, in this case, HA will not update the matrix A, so we skip this loop. 
		double norm_test = sqrt(cblas_ddot(m-k, sub_Av, 1, sub_Av, 1));
		if(sub_Av[0] == norm_test)
			continue;
		
		//hA is the reflector vector got from house function.
		double *hA = (double *) malloc((m-k)*sizeof(double));
		hA = house(sub_Av, (m-k), 1);
		
		//Copy the data in hA to vector v.
		for(i=0;i<m-k;i++)
		{
			v[i+k] = hA[i];
		}
		
		//Get the subset of A.
		double *sub_AM = (double *) malloc((m-k)*(n-k)*sizeof(double));
		for(row=0;row<m-k;row++)
		{
			for(col=0;col<n-k;col++)
			{
				sub_AM[row*(n-k)+col] = A[(row+k)*n+(col+k)];
			}
		}

		//Use subset of A and current reflector hA in applyH function to get the updated subset of A.
		double *aHAv = (double *) malloc((m-k)*(n-k)*sizeof(double));
		aHAv = applyH(sub_AM, m-k, n-k, hA, m-k, 1);

		//Copy the updated subset of A to A to update A.
		for(row=0;row<m-k;row++)
		{
			for(col=0;col<n-k;col++)
			{
				A[(row+k)*n+(col+k)] = aHAv[row*(n-k)+col];
			}
		}
		
		//Store the Householder vector in the strictly lower triangular part of A.
		for(i=k;i<m-1;i++)
			A[(i+1)*n+k] = v[i+1];
		
	}
	
	return A;
}
