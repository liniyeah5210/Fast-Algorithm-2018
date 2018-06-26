#include <stdio.h>
#include <stdlib.h> 
#include <math.h>

int Fast_Fourier_Transform(double *y_re, double *y_im, double *x_re, double *x_im, int N);

int main()
{
	int i, N;
	double *y_re, *y_im, *x_re, *x_im;
  	scanf ("%d",&N);	
	y_re = (double *) malloc( N * sizeof(double));
	y_im = (double *) malloc( N * sizeof(double));
	x_re = (double *) malloc( N * sizeof(double));
	x_im = (double *) malloc( N * sizeof(double));
	
	for(i=0;i<N;++i)
	{
		x_re[i] = i+1;
		x_im[i] = 0.0;
	}
	Fast_Fourier_Transform(y_re, y_im, x_re, x_im, N);
	for(i=0;i<N;++i)
	{
		printf("%f + %f i\n", y_re[i], y_im[i]);
	}
	
	 
}
int Fast_Fourier_Transform(double *y_re, double *y_im, double *x_re, double *x_im, int N)
{
		
	double w5re1, w5im1, w5re2, w5im2, w5re3, w5im3, w5re4, w5im4;
	w5re1 = cos(1*2*M_PI/5);
	w5im1 = -sin(1*2*M_PI/5);
	w5re2 = cos(2*2*M_PI/5); 
	w5im2 = -sin(2*2*M_PI/5);
	w5re3 = cos(3*2*M_PI/5); 
	w5im3 = -sin(3*2*M_PI/5); 
	w5re4 = cos(4*2*M_PI/5);
	w5im4 = -sin(4*2*M_PI/5);

	if(N==5) 
	{
		y_re[0] = x_re[0] + x_re[1] + x_re[2] + x_re[3] + x_re[4];
		y_im[0] = x_im[0] + x_im[1] + x_im[2] + x_im[3] + x_im[4];
		y_re[1] = x_re[0] + w3re1*x_re[1]-w3im1*x_im[1]+w3re2*x_re[2]-w3im2*x_im[2]; 
		y_im[1] = x_im[0] + w3re1*x_im[1]+w3im1*x_re[1]+w3re2*x_im[2]+w3im2*x_re[2];
		y_re[2] = x_re[0] + w3re2*x_re[1]-w3im2*x_im[1]+w3re1*x_re[2]-w3im1*x_im[2]; 
		y_im[2] = x_im[0] + w3re2*x_im[1]+w3im2*x_re[1]+w3re1*x_im[2]+w3im1*x_re[2];
	} else 
	{
		int k;
		double *y_0_re, *y_0_im, *y_1_re, *y_1_im, *y_2_re, *y_2_im;
		double *x_0_re, *x_0_im, *x_1_re, *x_1_im, *x_2_re, *x_2_im;
		double w_re, w_im, w2_re, w2_im, w_N_re, w_N_im, a1, a2, b1, b2, temp;
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
			x_0_re[k] = x_re[3*k];
			x_0_im[k] = x_im[3*k];
			x_1_re[k]  = x_re[3*k+1];
			x_1_im[k]  = x_im[3*k+1];
			x_2_re[k]  = x_re[3*k+2];
			x_2_im[k]  = x_im[3*k+2];
		}
		Fast_Fourier_Transform(y_0_re, y_0_im, x_0_re, x_0_im, N/3);
		Fast_Fourier_Transform(y_1_re, y_1_im, x_1_re, x_1_im, N/3);
		Fast_Fourier_Transform(y_2_re, y_2_im, x_2_re, x_2_im, N/3);
		

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

			y_re[k]     = y_0_re[k] + a1 + a2;
			y_im[k]     = y_0_im[k] + b1 + b2;
			y_re[k+N/3] = y_0_re[k] + a1*w3re1 - b1*w3im1 + a2*w3re2 - b2*w3im2;
			y_im[k+N/3] = y_0_im[k] + a1*w3im1 + b1*w3re1 + a2*w3im2 + b2*w3re2;  
			y_re[k+2*N/3] = y_0_re[k] + a1*w3re2 - b1*w3im2 + a2*w3re1 - b2*w3im1;  
			y_im[k+2*N/3] = y_0_im[k] + a1*w3im2 + b1*w3re2 + a2*w3im1 + b2*w3re1;    
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


	}
}
