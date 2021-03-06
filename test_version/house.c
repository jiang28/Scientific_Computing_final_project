//Compute the Householder vector v such that Hx = ||x||*e_1, where e_1 is the first
//unit vector and H = I - (2/v'v)vv'. The vector v is normalized so that v(1) =
//1.0. If ||x|| is already small, then an informational mesage is barfed out and
//a zero vector is returned..

#include <stdio.h>
#include <stdlib.h>
#include "cblas.h"
#include "math.h"

//Because there is no difference between row vector and column vector in C.
//Here I didn't follow the matlab code to receive a row vector in function house, instead, 
//I received the column vector sub_Av directly without change the receiving order of m and n.
double *house(double *x, int m, int n)	//m is the size of x, and n should be 1.
{

	//Set the parameter eps with a small enough number.
	double eps = 1E-6;
	int i;
        int inc = 1;
        double normx;

	double SMALL = sqrt(eps);

	if( m<1 || n!=1 )
	{
		printf("Size is wrong, Please reset the input vector x.");
	}
        
	normx = sqrt(cblas_ddot(m, x, inc, x, inc));
        printf("normx=%f \n" ,normx);
	if(normx < SMALL)
	{
		printf("Norm of x is effectively zero, so x is the zero vector already.\n");
		double *v = (double *) malloc(m*sizeof(double));
		for(i=0;i<m;i++)
			v[i] = 0;
		return v;
	}

	int sign;
	if(x[0]<0){
		sign = -1;}
	else{
		sign = 1;}
	
	//v = x - sign(x[0])*normx*e
	double *v = (double *) malloc(m*sizeof(double));

	//e is the first unit vector.
	double *e = (double *) malloc(m*sizeof(double));
	for(i=0;i<m;i++)
		e[i] = 0.0;
	e[0] = sign * normx;

	for(i=0;i<m;i++)
		v[i] = x[i] - e[i];
		
	//To normalize the first element to 1
	double first = v[0];
	if(first!=0)
	{
		for(i=0;i<m;i++)
			v[i] = v[i] / first;
	}
   
	return v;
}
