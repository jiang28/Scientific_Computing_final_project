#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include "cblas.h"
#include "applyH.h"
#include "applyQ.h"
#include "house.h"
#include "qrhouse.h"
#include <math.h>

char *input = "A_7_3";
char *output = "OUTPUT";

int main()
{
	int m;
	int n;
	int total;
	double *A;
	int row;
	int col;
	int i;
	
	//Read the data from file.
	FILE *fin;
	char line_buffer[BUFSIZ]; /* BUFSIZ is defined if you include stdio.h */
	char line_number;
	char * pEnd;
	fin = fopen(input, "r");
	long int li1;
	long int li2;
	double li3;
	line_number = 0;
	while (fgets(line_buffer, sizeof(line_buffer), fin))
	{
        	++line_number;
        	li1 = strtol (line_buffer,&pEnd,10);
        	li2 = strtol (pEnd,&pEnd,10);
        	li3 = strtod (pEnd,NULL);
        	int i1 = (int)li1;
        	int i2 = (int)li2;
		if(line_number == 1)
		{
			int i3 = (int)li3;
			m = i1;
			n = i2;
			total = i3;
			A = (double *) malloc(m*n*sizeof(double));
		}
		else
		{
			A[(i1-1)*n+(i2-1)] = li3;
		}
	}
	printf("Print matrix A read from file.\n");
	for(row=0;row<m;row++)
	{
		for(col=0;col<n;col++)
		{
			printf("A[%d][%d] = %lf\t", row, col, A[row*n+col]);
		}
		printf("\n");
	}


	struct timeval t_start;
	struct timeval t_end;
	gettimeofday(&t_start,NULL);
	printf("Start time is %d.\n",t_start.tv_usec);
	
	qrhouse(A, m, n);
	
	gettimeofday(&t_end,NULL);
	printf("End time is %d.\n",t_end.tv_usec);
	long timeA = t_end.tv_usec - t_start.tv_usec;
	printf("Total time is %dmicroseconds.\n",timeA);

	printf("Print matrix A after qrhouse function.\n");
	for(row=0;row<m;row++)
	{
		for(col=0;col<n;col++)
		{
			printf("A[%d][%d] = %lf\t", row, col, A[row*n+col]);
		}
		printf("\n");
	}
	
	
	//Solve b = applyQ(A, b)
	
	//Initialize b as an m size arithmetic sequence with 1, 2, 3,... m. 
	double *b_init = (double *) malloc(m*sizeof(double));
	for(i=0;i<m;i++)
		b_init[i] = 1.0/(i+1);
	
	double *b = (double *) malloc(m*sizeof(double));

	b = applyQ(A, m, n, b_init, m, 1);
	
	printf("Let's test b.\n");
	for(i=0;i<m;i++)
		printf("b[%d] = %lf\t", i, b[i]);
	printf("\n");
	
	
	//Solve x = triu(A(1:n,1:n))\b(1:n);
	
	//b_sub is the first n elements of b.
	double *b_sub = (double *) malloc(n*sizeof(double));
	for(i=0;i<n;i++)
		b_sub[i] = b[i];
	
	//A_triu is the upper triangular matrix of A.
	double *A_triu = (double *) malloc(n*n*sizeof(double));
	for(row=0;row<n;row++)
	{
		for(col=0;col<n;col++)
		{
			if(row<=col)
				A_triu[row*n+col] = A[row*n+col];
			else
				A_triu[row*n+col] = 0;
		}
	}
	
	//The result of A_triu\b_sub will be stored in b_sub.
	cblas_dtrsv(CblasRowMajor, CblasUpper, CblasNoTrans, CblasNonUnit, n, A_triu, n, b_sub, 1);
	
	printf("The final result x of the linear least squares problem.\n");
	for(i=0;i<n;i++)
	{
		printf("x[%d] = %lf\t", i, b_sub[i]);
	}
	printf("\n");
	

	//Print result to file.
	char cm[2] = "m";
	char cn[2] = "n";
	char ctime[20] = "time(microsecond)";
	char cflop_t[20] = "flops/microsecond";
	double vflop_t = (5*m*n*1.0)/timeA;
	printf("The flops/microsecond is %f.\n",vflop_t);

	FILE *fout;
	fout = fopen(output, "w");

	fprintf(fout, cm);
	fprintf(fout, "\t");
	fprintf(fout, cn);
	fprintf(fout, "\t");
	fprintf(fout, ctime);
	fprintf(fout, "\t");
	fprintf(fout, cflop_t);
	fprintf(fout, "\n");

	fprintf(fout, "%d", m);
	fprintf(fout, "\t");
	fprintf(fout, "%d", n);
	fprintf(fout, "\t");
	fprintf(fout, "%d", timeA);
	fprintf(fout, "\t\t\t");
	fprintf(fout, "%f", vflop_t);
	fprintf(fout, "\n");

	return 0;
}
