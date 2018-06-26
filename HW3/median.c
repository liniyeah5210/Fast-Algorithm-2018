#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <math.h>
#include <time.h>

#define DEBUG  1 //output all elements of array x and y
#define DEBUG2 0 //output the process in the selecting process

/*    HW3 find the median in the unsorted array
葉宜霖 C24031180


What's the median?
For a list with an odd number of elements the median is the (N/2)th smallest element.
For a list with an even number of elements the median is the average of the (N/2)th and (N/2-1)th elements.

*/

/*
The function "quickselect1" is the recursive function for finding the (N/2)th smallest element.
The function "quickselect2" is the recursive function for finding the (N/2-1)th smallest element after the (N/2)th smallest element is found by the function "quickselect1". So the function "quickselect2" is called in the implementation of the function "quickselect1".
*/


int quickselect1(int *x, int left, int right, int Ntot);
int quickselect2(int *x, int left, int right, int Ntot);

int main()
{
	clock_t t1, t2;				// variables for computing clocks 
	int *x, *y, s, p;
	double T1;
	int i, j, N;

		srand( time(NULL) );
		printf ("Enter the length of the array: ");
  		scanf ("%d",&N);	


		x = (int *) malloc( N * sizeof(int) );
		y = (int *) malloc( N * sizeof(int) );

		for(i=0;i<N;++i)
		{
			y[i] = x[i] = rand() % N;
		}
		#if DEBUG        // if DEBUG == 1, then compile the following codes 
		for(i=0;i<N;++i)
		{
			printf("x[%d]=%d\n",i,x[i]);
		}
		printf("\n");
		#endif			// end of if block

			
		for(i=0;i<N;++i) y[i] = x[i];
		
		t1 = clock();
		quickselect1(y,0,N,N);
		t2 = clock();
		T1 = (t2-t1)/(double) CLOCKS_PER_SEC;
		
		#if DEBUG
		printf("\n");
		for(i=0;i<N;++i)
		{
			printf("y[%d]=%d\n",i,y[i]);
		}
		#endif

			


		printf("The length of the array is: %d\n",N);
	
		if (N%2==1)
		printf("The median element is: y[%d]= %d\n",N/2, y[N/2]);
		if (N%2==0)
		printf("The median element is: (y[%d]+y[%d])/2 = (%d+%d)/2\n",N/2, N/2-1, y[N/2], y[N/2-1]);

		printf("The computation time for finding the median element: %f\n", T1);


			
		free(x);
		free(y);
	 

	return 0;
}



int quickselect1(int *x, int left, int right, int Ntot)
{
	int i, j, k, pivot, pivot_loc, N = right-left; 
	int *y;
	
	if(left < right-1)
	{
		y = (int *) malloc(N*sizeof(int));
		pivot_loc = left+(rand() % N);
		pivot = x[pivot_loc];
		x[pivot_loc] = x[left];
		x[left] = pivot;
		i = 0; j = N-1;
		 
		for(k=1;k<N;++k) 
		{
			if(x[left+k] <= pivot) 
			{
				y[i++] = x[left+k];          //  i = i + 1; 
				//i = i + 1;                 //  y[i] = x[left+k]; --> y[++i] = x[left+k];
			}
			else
			{
				y[j--] = x[left+k];
				// j = j - 1;
			}
		}
		y[i] = pivot;
		#if DEBUG2
		printf("%d %d %d %d %d\n",left,i,j,pivot,N);
		for(k=left;k<left+N;++k)
		{
			printf("y[%d]=%d\n",k,y[k]);
		}
		#endif
		for(k=0;k<N;++k)
		{
			x[left+k] = y[k];
		}
		free(y);
		
		//If the array has odd number of elements, find x[Ntot/2] 
		if (Ntot%2==1) {
			if (left+i<Ntot/2) quickselect1(x,left+i+1,right,Ntot);
			if (left+i>Ntot/2) quickselect1(x,left,left+i,Ntot);	
			if (left+i==Ntot/2) return 1;	
		}
		

		//If the array has even number of elements, find x[Ntot/2] and x[Ntot/2-1]
		//quickselect1 is for finding x[Ntot/2] and quickselect2 is for finding x[Ntot/2-1]
		if (Ntot%2==0) {
			if (left+i<Ntot/2) quickselect1(x,left+i+1,right,Ntot);
			if (left+i>Ntot/2) quickselect1(x,left,left+i,Ntot);	
			if (left+i==Ntot/2) 
			{	
				if (left+i!=left) quickselect2(x,left,left+i,Ntot);
				if (left+i==left) quickselect2(x,0,left+i,Ntot);
			}

		}



	}
	else
	{
		return 1;
	}
}


int quickselect2(int *x, int left, int right, int Ntot)
{
	int i, j, k, pivot, pivot_loc, N = right-left; 
	int *y;
	
	if(left < right-1)
	{
		y = (int *) malloc(N*sizeof(int));
		pivot_loc = left+(rand() % N);
		pivot = x[pivot_loc];
		x[pivot_loc] = x[left];
		x[left] = pivot;
		i = 0; j = N-1;
		 
		for(k=1;k<N;++k) 
		{
			if(x[left+k] <= pivot) 
			{
				y[i++] = x[left+k];          //  i = i + 1; 
				//i = i + 1;                 //  y[i] = x[left+k]; --> y[++i] = x[left+k];
			}
			else
			{
				y[j--] = x[left+k];
				// j = j - 1;
			}
		}
		y[i] = pivot;
		#if DEBUG2
		printf("%d %d %d %d %d\n",left,i,j,pivot,N);
		for(k=left;k<left+N;++k)
		{
			printf("y[%d]=%d\n",k,y[k]);
		}
		#endif
		for(k=0;k<N;++k)
		{
			x[left+k] = y[k];
		}
		free(y);
		
		

			if (left+i<Ntot/2-1) quickselect2(x,left+i+1,right,Ntot);
			if (left+i>Ntot/2-1) quickselect2(x,left,left+i,Ntot);
			if (left+i==Ntot/2) return 1;



	}



	
	else
	{
		return 1;
	}
}


