#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

/*
Fast Algorithm 2018
HW1
I-Lin Yeh (C24031180)
*/


int main()
{

	clock_t t1, t2;				// variables for computing clocks
	double a=1.234, b=2.456, c;
	double T1, T2;
	int i, j, k, N=100000000;

	// 1. compute the time for addition/subtraction
	t1 = clock();
	for(i=0;i<N;++i)
	{
		a = a+b;
		b = a-b;
		a = a-b;
	}
	t2 = clock();
	T1 = (t2-t1)/(double) CLOCKS_PER_SEC;
	printf("(+,-) time x 3 + 1 loop:%e\n",T1);
	printf("(a,b)=%f %f\n", a,b);
	
	t1 = clock();
	for(i=0;i<N;++i)
	{
		a = a+b;
		b = a-b;
		a = a-b;
		a = a+b;
		b = a-b;
		a = a-b;
	}
	t2 = clock();
	T2 = (t2-t1)/(double) CLOCKS_PER_SEC;
	printf("(+,-) time x 6 + 1 loop:%e\n",T2);
	printf("(a,b)=%f %f\n", a,b);
	printf("Real (+,-) time: %e\n",(T2-T1)/(3.0*N));
	printf("\n\n");


	// 2. compute the time for multiplication/division
	t1 = clock();
	for(i=0;i<N;++i)
	{
		a = a*b;
		b = a/b;
		a = a/b;
	}
	t2 = clock();
	T1 = (t2-t1)/(double) CLOCKS_PER_SEC;
	printf("(*,/) time x 3 + 1 loop:%e\n",T1);
	printf("(a,b)=%f %f\n", a,b);

	t1 = clock();
	for(i=0;i<N;++i)
	{
		a = a*b;
		b = a/b;
		a = a/b;
		a = a*b;
		b = a/b;
		a = a/b;
	}
	t2 = clock();
	T2 = (t2-t1)/(double) CLOCKS_PER_SEC;
	printf("(*,/) time x 6 + 1 loop:%e\n",T2);
	printf("(a,b)=%f %f\n", a,b);
	printf("Real (*,/) time:%e\n",(T2-T1)/(3.0*N));
	printf("\n\n");
	


	// 3. compute the time for sin
	t1 = clock();
	for(i=0;i<N;++i)
	{
		a = sin(a);
	}
	t2 = clock();
	T1 = (t2-t1)/(double) CLOCKS_PER_SEC;
	printf("(sin) time + 1 loop:%e\n",T1);
	printf("(a,b)=%f %f\n", a,b);

	t1 = clock();
	for(i=0;i<N;++i)
	{
		a = sin(b);
		b = sin(a);
	}
	t2 = clock();
	T2 = (t2-t1)/(double) CLOCKS_PER_SEC;
	printf("(sin) time x 2 + 1 loop:%e\n",T2);
	printf("(a,b)=%f %f\n", a,b);
	printf("Real (sin) time: %e\n",(T2-T1)/N);



	return 0;
}
