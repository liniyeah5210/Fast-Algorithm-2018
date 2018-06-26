#include <stdio.h>
#include <stdlib.h>
#include <math.h> 
#include <time.h>

/*
I use DST-I from wiki.
*/
int FFT(double *x_r, double *x_i, double *y_r, double *y_i, int N);
int DST(double *x_r, double *y_r, int N, int dir);

int main()
{
	int i, N, print, dir;
	double *x_r, *x_i, *x4_r, *x4_i, *y_r, *y_i; 
	clock_t t1, t2;


	// N must be 2^k-1	
	printf("N must be 2^k-1\n");
	printf("N=");
	scanf("%d", &N);
	printf("N=%d\n", N);
//	printf("dir=");
//	scanf("%d", &dir);
	printf("print=");
	scanf("%d", &print);
	
	x_r = (double *) malloc(N * sizeof(double));
	y_r = (double *) malloc(N * sizeof(double));

        for(i=0;i<N;++i)
	{
		x_r[i] = i; 
	}


/*
	for(i = 0; i < N; i ++)
	  {
	    printf("%2d> ",i);
	    scanf("%lf", &x_r[i]);
	  } 
*/

	DST(x_r, y_r, N, 1);
	DST(y_r, y_r, N, -1);





	if (print==1)	
	for(i=0;i<N;++i)
	{
		printf("y_r[%d] : %f \n", i, y_r[i]);
	}
	
	return 0;
}


int DST(double *x_r, double *y_r, int N, int dir)
{

	int i;
	double *x2_r,*x2_i,*y2_r,*y2_i; 
	x2_r = (double *) malloc(2*(N+1) * sizeof(double));
	x2_i = (double *) malloc(2*(N+1) * sizeof(double));
	y2_r = (double *) malloc(2*(N+1) * sizeof(double));
	y2_i = (double *) malloc(2*(N+1) * sizeof(double));
	

	if(dir!=1 && dir!=-1)
	{
	printf("dir must be +1 or -1\n");
	return 1;
	}

	for(i=0;i<2*(N+1);++i)
	{
		x2_r[i] = 0;
		x2_i[i] = 0;
		y2_r[i] = 0;
		y2_i[i] = 0;
	}
	
        for(i=0;i<N;++i)
	{
		x2_r[i+1] = x_r[i]; 
		x2_r[2*(N+1)-i-1] = -x_r[i];
	}
	

	FFT(x2_r, x2_i, y2_r, y2_i, 2*(N+1));



	if(dir==1)
	for(i=0;i<N;++i)
	{
		y_r[i] = -0.5*y2_i[i+1];
        }


	if(dir==-1)
	for(i=0;i<N;++i)
	{
		y_r[i] = -0.5*y2_i[i+1]*2/(N+1);
	}


}




int FFT(double *x_r, double *x_i, double *y_r, double *y_i, int N)
{
	int k;
	double w_re, w_im, w_N_re, w_N_im, temp, u1_r, u1_i, u2_r, u2_i, a, b, c, d;
	double *z_r, *z_i, *u_r, *u_i;

	if(N == 2)
	{
	y_r[0] = x_r[0] + x_r[1];
	y_i[0] = x_i[0] + x_i[1];
	y_r[1] = x_r[0] - x_r[1]; 
	y_i[1] = x_i[0] - x_i[1];
	} 


	else
	{
	double *z_r, *z_i, *u_r, *u_i;
	z_r = (double *) malloc(N*sizeof(double)); // z_r:0~N/2-1: even, z_r:N/2~N-1: odd
	if(z_r==NULL) { printf("no memory!!\n"); return 0; }
	z_i = (double *) malloc(N*sizeof(double)); // z_i:0~N/2-1: even, z_i:N/2~N-1: odd 
	if(z_i==NULL) { printf("no memory!!\n"); return 0; }
	u_r = (double *) malloc(N*sizeof(double)); // z_r:0~N/2-1: even, z_r:N/2~N-1: odd
	if(u_r==NULL) { printf("no memory!!\n"); return 0; }
	u_i = (double *) malloc(N*sizeof(double)); // z_i:0~N/2-1: even, z_i:N/2~N-1: odd 
	if(u_i==NULL) { printf("no memory!!\n"); return 0; }
	for(k=0;k<N/2;++k)
	{
		z_r[k] = x_r[2*k];
		z_i[k] = x_i[2*k];
		z_r[N/2+k]  = x_r[2*k+1];
		z_i[N/2+k]  = x_i[2*k+1];
	}		
	FFT(z_r, z_i, u_r, u_i, N/2);
	FFT(z_r+N/2, z_i+N/2, u_r+N/2, u_i+N/2, N/2);
	w_N_re =  cos(2.0*M_PI/N);
	w_N_im = -sin(2.0*M_PI/N);
	w_re   = 1.0;
	w_im   = 0.0; 
	for(k=0;k<N/2;++k)
	{
		a = w_re*u_r[N/2+k] - w_im*u_i[N/2+k];
		b = w_re*u_i[N/2+k] + w_im*u_r[N/2+k];
		y_r[k]     = u_r[k] + a;
		y_i[k]     = u_i[k] + b;
		y_r[N/2+k] = u_r[k] - a;
		y_i[N/2+k] = u_i[k] - b;
		temp = w_re;
		w_re = w_re*w_N_re - w_im*w_N_im;
		w_im = temp*w_N_im + w_im*w_N_re;
	}
	free(u_r); free(u_i); free(z_r); free(z_i);
	}	
	return 0;
} 
