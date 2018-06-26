#include <stdio.h>
#include <stdlib.h>
#include <math.h> 
#include <time.h>

int main()
{
	int i, N, N2, K, *p;
	double *x_r, *x_i, *y_r, *y_i; // y = fft(x);
	clock_t t1, t2;
	
	printf("N=");
	scanf("%d", &N);
	printf("N=%d\n", N);
	
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
//	system("pause");
//	for(i=0;i<N;++i)
//	{
//		printf("y[%d] : %f + %f i\n", i, y_r[i], y_i[i]);
//	}

	return 0;
} 
int FFT_general(double *x_r, double *x_i, double *y_r, double *y_i, int N)
{
	int k, j, p;
	double w_re, w_im, w_N_re, w_N_im, temp, u1_r, u1_i, u2_r, u2_i, a, b, c, d;
	double a_r[5], a_i[5], b_r[5], b_i[5]; 

	if(N == 2 || N == 3 || N == 5)
	{
		DFT(y_r, y_i, x_r, x_i, N);
	} else 
	{
		double *z_r, *z_i, *u_r, *u_i;
		z_r = (double *) malloc(N*sizeof(double)); // z_r:0~N/2-1: even, z_r:N/2~N-1: odd
		z_i = (double *) malloc(N*sizeof(double)); // z_i:0~N/2-1: even, z_i:N/2~N-1: odd 
		u_r = (double *) malloc(N*sizeof(double)); // z_r:0~N/2-1: even, z_r:N/2~N-1: odd
		u_i = (double *) malloc(N*sizeof(double)); // z_i:0~N/2-1: even, z_i:N/2~N-1: odd 
		
		if(N % 2 == 0) p = 2;
		else if(N % 3 == 0) p = 3;
		else if(N % 5 == 0) p = 5;
		else {
			p = N;
			DFT(y_r, y_i, x_r, x_i, N);
			return 0;
		}

		for(k=0;k<N/p;++k)
		{
			for(j=0;j<p;++j)
			{
				z_r[j*N/p+k] = x_r[p*k+j];
				z_i[j*N/p+k] = x_i[p*k+j];
			}
		}
		for(j=0;j<p;++j)
		{
			FFT_general(z_r+j*N/p,z_i+j*N/p,u_r+j*N/p,u_i+j*N/p,N/p);
		}
		w_N_re =  cos(2.0*M_PI/N);
		w_N_im = -sin(2.0*M_PI/N);
		w_re   = 1.0;
		w_im   = 0.0; 
		for(k=0;k<N/p;++k)
		{
			for(j=0;j<p;++j) 
			{
				a_r[j] = u_r[j*N/p+k];
				a_i[j] = u_i[j*N/p+k];
			}
			DFT_wk(a_r,a_i,w_re,w_im,p);
			DFT(b_r,b_i,a_r,a_i,p);
			for(j=0;j<p;++j) 
			{
				y_r[j*N/p+k] = b_r[j];
				y_i[j*N/p+k] = b_i[j];
			}
			temp = w_re;
			w_re = w_re*w_N_re - w_im*w_N_im;
			w_im = temp*w_N_im + w_im*w_N_re;
		}
		
		/*
		if(N % 2 == 0)
		{

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
				for(j=0;j<2;++j) 
				{
					a_r[j] = u_r[j*N/2+k];
					a_i[j] = u_i[j*N/2+k];
				}
				DFT_wk(a_r,a_i,w_re,w_im,2);
				DFT(b_r,b_i,a_r,a_i,2);
				for(j=0;j<2;++j) 
				{
					y_r[j*N/2+k] = b_r[j];
					y_i[j*N/2+k] = b_i[j];
				}				
				temp = w_re;
				w_re = w_re*w_N_re - w_im*w_N_im;
				w_im = temp*w_N_im + w_im*w_N_re;
			}
		} else if(N % 3 == 0)
		{
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
				for(j=0;j<3;++j) 
				{
					a_r[j] = u_r[j*N/3+k];
					a_i[j] = u_i[j*N/3+k];
				}
				DFT_wk(a_r,a_i,w_re,w_im,3);
				DFT(b_r,b_i,a_r,a_i,3);
				for(j=0;j<3;++j) 
				{
					y_r[j*N/3+k] = b_r[j];
					y_i[j*N/3+k] = b_i[j];
				}
				temp = w_re;
				w_re = w_re*w_N_re - w_im*w_N_im;
				w_im = temp*w_N_im + w_im*w_N_re;
			}
		} else if(N % 5 == 0)
		{
			for(k=0;k<N/5;++k)
			{
				for(j=0;j<5;++j)
				{
					z_r[j*N/5+k] = x_r[5*k+j];
					z_i[j*N/5+k] = x_i[5*k+j];
				}
			}
			for(j=0;j<5;++j)
			{
				FFT_general(z_r+j*N/5,z_i+j*N/5,u_r+j*N/5,u_i+j*N/5,N/5);
			}
			w_N_re =  cos(2.0*M_PI/N);
			w_N_im = -sin(2.0*M_PI/N);
			w_re   = 1.0;
			w_im   = 0.0; 
			for(k=0;k<N/5;++k)
			{
				for(j=0;j<5;++j) 
				{
					a_r[j] = u_r[j*N/5+k];
					a_i[j] = u_i[j*N/5+k];
				}
				DFT_wk(a_r,a_i,w_re,w_im,5);
				DFT(b_r,b_i,a_r,a_i,5);
				for(j=0;j<5;++j) 
				{
					y_r[j*N/5+k] = b_r[j];
					y_i[j*N/5+k] = b_i[j];
				}
				temp = w_re;
				w_re = w_re*w_N_re - w_im*w_N_im;
				w_im = temp*w_N_im + w_im*w_N_re;
			}
		}
		*/
		free(u_r); free(u_i); free(z_r); free(z_i);
		
	}	
	return 0;
}
int DFT(double *y_r, double *y_i, double *x_r, double *x_i, int N)
{
	int k, n;
	double a, c, s, ca, sa, t;
	
	for(k=0;k<N;++k)
	{
		y_r[k] = 0.0;
		y_i[k] = 0.0;
		a  = 2*M_PI*k/N;
		ca = cos(a);
		sa = sin(a);
		c  = 1.0;
		s  = 0.0;
		for(n=0;n<N;++n)
		{
			y_r[k] += x_r[n]*c + x_i[n]*s;
			y_i[k] += x_i[n]*c - x_r[n]*s;
			t =  c;
			c =  c*ca - s*sa;
			s =  s*ca + t*sa; 
		}
	}
	return 0;
}

int DFT_wk(double *x_r, double *x_i, double w_r, double w_i, int N)
{
	int k;
	double c, s, t;

	c  = 1.0;
	s  = 0.0;	
	for(k=1;k<N;++k)
	{
		t = c;
		c = c*w_r-s*w_i;
		s = t*w_i+s*w_r;
		t = x_r[k];
		x_r[k] = x_r[k]*c-x_i[k]*s;
		x_i[k] = t     *s+x_i[k]*c;
	}
	return 0;
}
