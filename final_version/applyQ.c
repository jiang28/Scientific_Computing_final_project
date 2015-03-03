// Multiply an orthogonal matrix times a vector b, where the orthogonal matrix
// is implicitly defined as a product of Householder reflectors, with A 
// containing the Householder vector.

#include <stdio.h>
#include <stdlib.h>
#include "cblas.h"
#include "applyH.h"
#include <math.h>

double *applyQ(double *A, int mA, int nA, double *b, int mb, int nb)
{
	if(mb!=mA || nb>1)
	{
		printf("Size is wrong, Please reset the input matrix A and vector b.");
		double *p = NULL;
		A = p;
		return b;
	}

	int m = mA;
	int n = nA;
	int i;

	int k;
	for(k=0;k<n;k++)
	{
		//The array v is used to retrieve back the Householder reflector stored in A.
		double *v = (double *) malloc((m-k)*sizeof(double));
		v[0] = 1;
		//This tag is used to test whether this reflector vector is all-zero vector.
		int tag = 0;
		for(i=1;i<m-k;i++)
		{
			v[i] = A[(i+k)*n+k];
			if(v[i]!=0)
				tag = 1;
		}
		if(tag==0)
			continue;
		//Because if the sub_Av is a vector that only has non-zero element in the first position, the reflector will be all zeros.
		//Thus the reflector should be retrieved as all-zero, instead of 1 in the first position. In this case, the H will be
		//an identity matrix, which will not change the b vector in applyH function below, so we skip this loop.

		double *aHbv = (double *) malloc((m-k)*sizeof(double));
		double *sub_b = (double *) malloc((m-k)*sizeof(double));
		for(i=0;i<m-k;i++)
			sub_b[i] = b[i+k];
		
		aHbv = applyH(sub_b, (m-k), 1, v, (m-k), 1);
		
		for(i=0;i<m-k;i++)
			b[i+k] = aHbv[i];
	}
	return b;
}
