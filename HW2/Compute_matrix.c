#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <omp.h>

int main()
{
	clock_t t1, t2;				 
	double *x, *b, *c,T1, T2, T3;		
	double **A;
	int i, j, k, L, M, N=1000;		
		
		A = (double **) malloc( N * sizeof(double *));  
		A[0] = (double *) malloc(N*N*sizeof(double));

		#pragma omp parallel for
		for(i=1;i<N;++i) 
		{
			A[i] = A[0] + i*N;
		}
		// A[1] = A[0] + N, A[2] = A[0]+2N, A[3] = A[0]+3N  

		x = (double *) malloc( N * sizeof(double));
		b = (double *) malloc( N * sizeof(double));
		c = (double *) malloc( N * sizeof(double));
		
		//Generate matrix parallelly 
		M = N/4;
		#pragma omp parallel num_threads(4) private(i,j,k,L)
		{
			k = omp_get_thread_num();
			srand(time(NULL)>>k);			 
			for(i=k*M;i<(k+1)*M;++i)
			{
				for(j=0;j<N;++j)
				{
					A[i][j] = rand() % 10;
				}
				x[i] = rand();
			}
		}	


		//Compute matrix-vector multiplication without parallelization
		double t;
		t1 = clock();
		for(i=0;i<N;++i) 
		{
			t = 0.0;
			for(j=0;j<N;++j)
			{
				t += A[i][j]*x[j];
			}
			b[i] = t;
		}
		t2 = clock();
		T1 = (t2-t1)/(double)CLOCKS_PER_SEC;

		
		//Compute matrix-vector multiplication with parallelization 
		t1 = clock();
		M = N/10;
		#pragma omp parallel num_threads(10) private(i,j,k,t)
		{
			k = omp_get_thread_num();
			for(i=k*M;i<(k+1)*M;++i) 
			{
				t = 0.0;
				for(j=0;j<N;++j)
				{
					t += A[i][j]*x[j];
				}
				c[i] = t;
			}
		}
		t2 = clock();
		T2 = (t2-t1)/(double) CLOCKS_PER_SEC;


		


		//Compare the answers and print the time
		for(i=0;i<N;++i)
		{
			if(fabs(b[i]-c[i])!=0)
			printf("Parallelization method 1: Wrong at i = %d, the fabs(b[i]-c[i]) = %f\n", i, fabs(b[i]-c[i]));
		}
		printf("Matrix time vector : %f\n",T1);
		printf("Matrix time vector (parallel) : %f\n",T2);
		free(b);
		free(x);
		free(A[0]);
		free(A);

	return 0;
} 
