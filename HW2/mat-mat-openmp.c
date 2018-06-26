#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>
//#include <omp.h>

#define  N      1000
#define  DEBUG  1
#define  EPS    1.0e-18

/* Please define here for the matrices. */
double  A[N][N];
double  B[N][N];
double  C[N][N];

void MyMatMat(double C[N][N], double A[N][N], double B[N][N], int n); 

int main(int argc, char* argv[]) {

     double  t1, t2, t_w;
     double  dc_inv, d_mflops;

     int     i, j;      
     int     iflag;

     int    myid, numprocs;
     int    ierr, rc;


     ierr = MPI_Init(&argc, &argv);
     ierr = MPI_Comm_rank(MPI_COMM_WORLD, &myid);
     ierr = MPI_Comm_size(MPI_COMM_WORLD, &numprocs);


     /* matrix generation --------------------------*/
     if (DEBUG == 1) {
       for(j=0; j<N; j++) {
         for(i=0; i<N; i++) {
           A[j][i] = 1.0;
           B[j][i] = 1.0;
           C[j][i] = 0.0;
         }
       }
     } else {
       srand(1);
       dc_inv = 1.0/(double)RAND_MAX;
       for(j=0; j<N; j++) {
         for(i=0; i<N; i++) {
           A[j][i] = rand()*dc_inv;
           B[j][i] = rand()*dc_inv;
           C[j][i] = 0.0;
         }
       }
     } /* end of matrix generation --------------------------*/

     /* Start of mat-vec routine ----------------------------*/
     t1 = MPI_Wtime();

     MyMatMat(C, A, B, N);

     t2 = MPI_Wtime();
     t_w =  t2 - t1; 
     /* End of mat-vec routine --------------------------- */

     printf("N  = %d \n",N);
     printf("Mat-Mat time  = %lf [sec.] \n",t_w);

     d_mflops = 2.0*(double)N*(double)N*(double)N/t_w;
     d_mflops = d_mflops * 1.0e-6;
     printf(" %lf [MFLOPS] \n", d_mflops);

     /* Verification routine ----------------- */
     if (DEBUG == 1) {
       iflag = 0;
       for(j=0; j<N; j++) { 
         for(i=0; i<N; i++) { 
           if (fabs(C[j][i] - (double)N) > EPS) {
             printf(" Error! in ( %d , %d ) th argument. \n",j, i);
             iflag = 1;
             break;
           } 
         }
       }
       if (iflag == 0) printf("OK! \n");
     }
     /* ------------------------------------- */

     rc = MPI_Finalize();

}

void MyMatMat(double C[N][N], double A[N][N], double B[N][N], int n) 
{
     int  i, j, k;
#pragma omp parallel for private (j, k)
     for(i=0; i<n; i++) {
       for(j=0; j<n; j++) {
         for(k=0; k<n; k++) {
           C[i][j] += A[i][k] * B[k][j]; 
         }
       }
     }
   
}

