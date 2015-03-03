#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include "cblas.h"
#include "applyH.h"
#include "applyQ.h"
#include "house.h"
#include "qrhouse.h"
#include <math.h>

//You can change the number of files and file list below.
int number_of_files = 12;
char *input[] = {"A_4_3", "A_7_3", "A_17_4", "A_35_6", "A_40_10", "A_40_30", "A_50_40", "A_80_10", "A_100_80", "A_250_200", "A_500_400", "A_1000_800"};
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
	
	//Print the first line to file.
	char cm[2] = "m";
	char cn[2] = "n";
	char ctime[20] = "time";
	char cflop_t[20] = "gigaFlops/s";

	FILE *fout;
	fout = fopen(output, "w");

	fprintf(fout, cm);
	fprintf(fout, "\t");
	fprintf(fout, cn);
	fprintf(fout, "\t");
	fprintf(fout, ctime);
	fprintf(fout, "\t\t\t");
	fprintf(fout, cflop_t);
	fprintf(fout, "\n");
	
	int k;
	for(k=0;k<number_of_files;k++)
	{
		printf("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %s ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n", input[k]);
		
		//Read the data from file.
		FILE *fin;
		char line_buffer[BUFSIZ]; /* BUFSIZ is defined if you include stdio.h */
		long int line_number;
		char * pEnd;
		fin = fopen(input[k], "r");
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

		struct timeval t_start;
		struct timeval t_end;
		gettimeofday(&t_start,NULL);
		printf("Start time is %d.%d\n", t_start.tv_sec, t_start.tv_usec);
	
		qrhouse(A, m, n);
	
		gettimeofday(&t_end,NULL);
		printf("End time is %d.%d\n", t_end.tv_sec, t_end.tv_usec);
		long timeA = (t_end.tv_sec - t_start.tv_sec) * 1E6 + (t_end.tv_usec - t_start.tv_usec);
		double time = ((double)(timeA))/1E6;
		printf("Total time is %10.17lf s.\n",time);
		
	
		//Solve b = applyQ(A, b)
	
		//Initialize b as an m size arithmetic sequence with 1, 2, 3,... m. 
		double *b_init = (double *) malloc(m*sizeof(double));
		for(i=0;i<m;i++)
			b_init[i] = 1.0/(i+1);
	
		double *b = (double *) malloc(m*sizeof(double));

		b = applyQ(A, m, n, b_init, m, 1);
		
	
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
			printf("x[%d] = %lf ", i, b_sub[i]);
		}
		printf("\n\n");
		
	
		double vflop_t = ((5*m*n*1.0)/1E9)/time;

		fprintf(fout, "%d", m);
		fprintf(fout, "\t");
		fprintf(fout, "%d", n);
		fprintf(fout, "\t");
		fprintf(fout, "%10.17lf", time);
		fprintf(fout, "\t");
		fprintf(fout, "%10.17lf", vflop_t);
		fprintf(fout, "\n");
	}
	
	return 0;
}
