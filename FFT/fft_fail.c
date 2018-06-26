#include <stdio.h>
#include <stdlib.h>
#include <math.h> 
#include <time.h>

int main()
{
	int i, N, N2, K,print, *p;
	double *x_r, *x_i, *y_r, *y_i; // y = fft(x);
	clock_t t1, t2;
	
	printf("N=");
	scanf("%d", &N);
	printf("N=%d\n", N);
	printf("print=");
	scanf("%d", &print);
	
	x_r = (double *) malloc(N * sizeof(double));
	x_i = (double *) malloc(N * sizeof(double));
	y_r = (double *) malloc(N * sizeof(double));
	y_i = (double *) malloc(N * sizeof(double));
	
	for(i=0;i<N;++i)
	{
		x_r[i] = i+1;
		x_i[i] = 0;
	}
	t1 = clock();
	FFT_general(x_r, x_i, y_r, y_i, N);
	t2 = clock();
	printf("time=%f\n",1.0*(t2-t1)/CLOCKS_PER_SEC);
	if (print==1)	
	for(i=0;i<N;++i)
	{
		printf("y[%d] : %f + %f i\n", i, y_r[i], y_i[i]);
	}
	
	return 0;
} 
int FFT_general(double *x_r, double *x_i, double *y_r, double *y_i, int N)
{
	int k;
	double w_re, w_im, w_N_re, w_N_im, temp, u1_r, u1_i, u2_r, u2_i, a, b, c, d;

	double w3re1, w3im1, w3re2, w3im2;
	w3re1 = cos(1*2*M_PI/3);
	w3im1 = -sin(1*2*M_PI/3);
	w3re2 = cos(2*2*M_PI/3); 
	w3im2 = -sin(2*2*M_PI/3);
	/*
	p = (int *) malloc(N * sizeof(int));	
	N2 = N;
	K  = 0;
	while(N2>1)
	{
		if(N2 % 2 == 0) { N2 /= 2; p[K] = 2; }
		else if(N2 % 3 == 0) { N2 /= 3; p[K] = 3; }
		else if(N2 % 5 == 0) { N2 /= 5; p[K] = 5; }
		K++;
	}
	printf("K=%d\n", K); 
	for(i=0;i<K;++i)
	{
		printf("%d ", p[i]);
	}
	*/
	if(N == 2)
	{
		y_r[0] = x_r[0] + x_r[1];
		y_i[0] = x_i[0] + x_i[1];
		y_r[1] = x_r[0] - x_r[1]; 
		y_i[1] = x_i[0] - x_i[1];
	} else if (N == 3)
	{
		// w_3    = cos(2pi/3)-i sin(2pi/3)=(-1/2 - i*sqrt(3)/2) = a+bi
		// w_3^2  = (-1/2 + i*sqrt(3)/2) = a-bi
		// y[0] = x[0]+x[1]+x[2]
		// y[1] = x[0]+x[1]*w_3+x[2]*w_3^2
		// y[2] = x[0]+x[1]*w_3^2 + x[2]*w_3
		y_r[0] = x_r[0] + x_r[1] + x_r[2];
		y_i[0] = x_i[0] + x_i[1] + x_i[2];
		y_r[1] = x_r[0] + w3re1*x_r[1]-w3im1*x_i[1]+w3re2*x_r[2]-w3im2*x_i[2]; 
		y_i[1] = x_i[0] + w3re1*x_i[1]+w3im1*x_r[1]+w3re2*x_i[2]+w3im2*x_r[2];
		y_r[2] = x_r[0] + w3re2*x_r[1]-w3im2*x_i[1]+w3re1*x_r[2]-w3im1*x_i[2]; 
		y_i[2] = x_i[0] + w3re2*x_i[1]+w3im2*x_r[1]+w3re1*x_i[2]+w3im1*x_r[2];
	} else if (N == 5)
	{
		y_r[0] = x_r[0] + x_r[1];
		y_i[0] = x_i[0] + x_i[1];
		y_r[1] = x_r[0] - x_r[1]; 
		y_i[1] = x_i[0] - x_i[1];
		y_r[2] = x_r[0] - x_r[1]; 
		y_i[2] = x_i[0] - x_i[1];	
		y_r[3] = x_r[0] - x_r[1]; 
		y_i[3] = x_i[0] - x_i[1];	
		y_r[4] = x_r[0] - x_r[1]; 
		y_i[4] = x_i[0] - x_i[1];	
	} else 
	{
		if(N % 2 == 0)
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
			FFT_general(z_r, z_i, u_r, u_i, N/2);
			FFT_general(z_r+N/2, z_i+N/2, u_r+N/2, u_i+N/2, N/2);
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
		} else if(N % 3 == 0)
		{
/*			double *z_r, *z_i, *u_r, *u_i;
			z_r = (double *) malloc(N*sizeof(double)); // z_r:0~N/2-1: even, z_r:N/2~N-1: odd
			z_i = (double *) malloc(N*sizeof(double)); // z_i:0~N/2-1: even, z_i:N/2~N-1: odd 
			u_r = (double *) malloc(N*sizeof(double)); // z_r:0~N/2-1: even, z_r:N/2~N-1: odd
			u_i = (double *) malloc(N*sizeof(double)); // z_i:0~N/2-1: even, z_i:N/2~N-1: odd 

			for(k=0;k<N/3;++k)
			{
				z_r[k] = x_r[3*k];
				z_i[k] = x_i[3*k];
				z_r[N/3+k]  = x_r[3*k+1];
				z_i[N/3+k]  = x_i[3*k+1];
				z_r[2*N/3+k]  = x_r[3*k+2];
				z_i[2*N/3+k]  = x_i[3*k+2];
			}
			FFT_general(z_r, z_i, u_r, u_i, N/3);
			FFT_general(z_r+N/3, z_i+N/3, u_r+N/3, u_i+N/3, N/3);
			FFT_general(z_r+2*N/3, z_i+2*N/3, u_r+2*N/3, u_i+2*N/3, N/3);
			w_N_re =  cos(2.0*M_PI/N);
			w_N_im = -sin(2.0*M_PI/N);
			w_re   = 1.0;
			w_im   = 0.0; 
			for(k=0;k<N/3;++k)
			{
				u1_r = w_re*u_r[N/3+k] - w_im*u_i[N/3+k];
				u1_i = w_re*u_i[N/3+k] + w_im*u_r[N/3+k];
				// (w_re + i w_im)^2 = w_re*w_re - w_im*w_im + i (2*w_re*w_im)
				u2_r = (w_re*w_re - w_im*w_im)*u_r[2*N/3+k] - (2*w_re*w_im)*u_i[2*N/3+k];
				u2_i = (w_re*w_re - w_im*w_im)*u_i[2*N/3+k] + (2*w_re*w_im)*u_r[2*N/3+k];
				a = -0.5; b = -sqrt(3)/2;
				y_r[k]       = u_r[k] + u1_r + u2_r;
				y_i[k]       = u_i[k] + u1_i + u2_i;
				y_r[N/3+k]   = u_r[k] + (u1_r*a-u1_i*b) + (u2_r*a+u2_i*b);
				y_i[N/3+k]   = u_i[k] + (u1_i*a+u1_r*b) + (u2_i*a-u2_r*b);
				y_r[2*N/3+k] = u_r[k] + (u1_r*a+u1_i*b) + (u2_r*a-u2_i*b);
				y_i[2*N/3+k] = u_i[k] + (u1_i*a-u1_r*b) + (u2_i*a+u2_r*b);
				temp = w_re;
				w_re = w_re*w_N_re - w_im*w_N_im;
				w_im = temp*w_N_im + w_im*w_N_re;
			}
			free(u_r); free(u_i); free(z_r); free(z_i);
		
*/

		double *y_0_re, *y_0_im, *y_1_re, *y_1_im, *y_2_re, *y_2_im;
		double *x_0_re, *x_0_im, *x_1_re, *x_1_im, *x_2_re, *x_2_im;
		double w_re, w_im, w2_re, w2_im, a1, a2, b1, b2;
		y_0_re = (double *) malloc( N/3 * sizeof(double));
		y_0_im = (double *) malloc( N/3 * sizeof(double));
		x_0_re = (double *) malloc( N/3 * sizeof(double));
		x_0_im = (double *) malloc( N/3 * sizeof(double));
		y_1_re = (double *) malloc( N/3 * sizeof(double));
		y_1_im = (double *) malloc( N/3 * sizeof(double));
		x_1_re = (double *) malloc( N/3 * sizeof(double));
		x_1_im = (double *) malloc( N/3 * sizeof(double));
		y_2_re = (double *) malloc( N/3 * sizeof(double));
		y_2_im = (double *) malloc( N/3 * sizeof(double));
		x_2_re = (double *) malloc( N/3 * sizeof(double));
		x_2_im = (double *) malloc( N/3 * sizeof(double));
		for(k=0;k<N/3;++k)
		{
			x_0_re[k] = x_r[3*k];
			x_0_im[k] = x_i[3*k];
			x_1_re[k]  = x_r[3*k+1];
			x_1_im[k]  = x_i[3*k+1];
			x_2_re[k]  = x_r[3*k+2];
			x_2_im[k]  = x_i[3*k+2];
		}
		FFT_general(y_0_re, y_0_im, x_0_re, x_0_im, N/3);
		FFT_general(y_1_re, y_1_im, x_1_re, x_1_im, N/3);
		FFT_general(y_2_re, y_2_im, x_2_re, x_2_im, N/3);
		

		w_N_re =  cos(2.0*M_PI/N);
		w_N_im = -sin(2.0*M_PI/N);
		w_re   = 1.0;
		w_im   = 0.0; 
		for(k=0;k<N/3;++k)
		{
			a1 = w_re*y_1_re[k] - w_im*y_1_im[k];
			b1 = w_re*y_1_im[k] + w_im*y_1_re[k];

			w2_re = w_re*w_re - w_im*w_im;
			w2_im = w_re*w_im + w_im*w_re;	

			a2 = w2_re*y_2_re[k] - w2_im*y_2_im[k];
			b2 = w2_re*y_2_im[k] + w2_im*y_2_re[k];

			y_r[k]     = y_0_re[k] + a1 + a2;
			y_i[k]     = y_0_im[k] + b1 + b2;
			y_r[k+N/3] = y_0_re[k] + a1*w3re1 - b1*w3im1 + a2*w3re2 - b2*w3im2;
			y_i[k+N/3] = y_0_im[k] + a1*w3im1 + b1*w3re1 + a2*w3im2 + b2*w3re2;  
			y_r[k+2*N/3] = y_0_re[k] + a1*w3re2 - b1*w3im2 + a2*w3re1 - b2*w3im1;  
			y_i[k+2*N/3] = y_0_im[k] + a1*w3im2 + b1*w3re2 + a2*w3im1 + b2*w3re1;    
			temp = w_re;
			w_re = w_re*w_N_re - w_im*w_N_im;
			w_im = temp*w_N_im + w_im*w_N_re;
		}
		free(y_0_re);
		free(x_0_re);
		free(y_0_im);
		free(x_0_im);
		free(y_1_re);
		free(y_1_im);
		free(x_1_re);
		free(x_1_im);
		free(y_2_re);
		free(y_2_im);
		free(x_2_re);
		free(x_2_im);









} else if(N % 5 == 0)
		{
			
		}
	}	
	//free(p);
	return 0;
} 
