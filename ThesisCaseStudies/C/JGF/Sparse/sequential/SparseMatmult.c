#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <time.h>
#include "SparseMatmult.h"


void run (int size, int validation){
        
    Sparse *sparse = malloc(sizeof(Sparse));
    sparse->size = size;
    
    if(!JGFinitialise(sparse))
    {
        printf("Memory problems at JGFinitialise! \n");
    }
    else
    {
        JGFkernel(sparse);

        if(validation)
        {
            JGFvalidate(sparse);
        }
    
        free_sparse(&sparse);
    }
}

void free_sparse(Sparse **sparse){
    
    free((*sparse)->col);
    free((*sparse)->row);
    free((*sparse)->val);
    free((*sparse)->x);
    free((*sparse)->y);
    free((*sparse));
    *sparse         = NULL;
}

int JGFinitialise(Sparse *sparse){

    int datasizes_M[]   = { 50000, 100000,  500000, 1000000, 1500000 };  
    int datasizes_N[]   = { 50000, 100000,  500000, 1000000, 1500000 };
    int datasizes_nz[]  = {250000, 500000, 2500000, 5000000, 7500000 };
   
    int size                = sparse->size;
    
    sparse->datasizes_M     = datasizes_M[size];
    sparse->datasizes_N     = datasizes_N[size];
    sparse->datasizes_nz    = datasizes_nz[size];
    
    sparse->SPARSE_NUM_ITER = 200;
    sparse->RANDOM_SEED     = 10101010;
    
    srand(sparse->RANDOM_SEED);
    
    sparse->x       = RandomVector(sparse->datasizes_N);
    sparse->y       = malloc(sizeof(double) * sparse->datasizes_M);
    sparse->val     = malloc(sizeof(double) * sparse->datasizes_nz);
    sparse->col     = malloc(sizeof(int)    * sparse->datasizes_nz);
    sparse->row     = malloc(sizeof(int)    * sparse->datasizes_nz);
    
    
    
    if(!(sparse->x  && sparse->y && 
            sparse->val && sparse->col && sparse->row))
    {
        return 0;
    }
    
    init_double_array(sparse->datasizes_M, sparse->y);
    
    for (int i = 0; i < sparse->datasizes_nz; i++) 
    {
        // generate random row index (0, M-1)
        sparse->row[i] = abs_(rand()) % sparse->datasizes_M;

        // generate random column index (0, N-1)
        sparse->col[i] = abs_(rand()) % sparse->datasizes_N;

        sparse->val[i] = (((double)rand()/(double)RAND_MAX));
    }
    return 1;
  }

double *RandomVector(int N){
    
    double *A = malloc(sizeof(double) * N);
    
    for (int i = 0; i < N; i++)
    {
        A[i] = (((double)rand()/(double)RAND_MAX)) * 1e-6;
    }

    return A;
}

int abs_            (int value)     {return (value > 0) ? value : -value;}
double abs_double   (double value)  {return (value > 0) ? value : -value;}

void JGFkernel(Sparse *sparse){
    
    sparse->ytotal = test(sparse->y, sparse->val, sparse->row,
                         sparse->col, sparse->x, sparse->SPARSE_NUM_ITER,
                         sparse->datasizes_nz);    
}

double test (   double y[], double val[], int row[],
                int col[], double x[], int NUM_ITERATIONS, 
                int var_length
            ){
    
    int nz = var_length;
    double ytotal = 0;
    
    const double start  = omp_get_wtime();
    
    for (int reps = 0; reps < NUM_ITERATIONS; reps++)
    {
        for (int i = 0; i < nz; i++)
        {
            y[ row[i] ] += x[ col[i] ] * val[i];
        }
    }
    
    const double end    = omp_get_wtime();
    printf("%f\n", (end - start));

    for (int i = 0; i < nz; i++) 
    {
        ytotal += y[ row[i] ];
    }
    return ytotal;
}

void JGFvalidate(Sparse *sparse){
    
    double refval[] = { 75.205209621081394, 150.325289707762209, 
                        750.707296942218136, 1499.080290091720144, 
                        2250.518217459347852};
    
    double dev = abs_double(sparse->ytotal - refval[sparse->size]);
    if (dev > 1.0e-12)
    {
        printf("Validation failed \n");
	printf("ytotal = %.15f %.15f %d \n", sparse->ytotal, dev, sparse->size);
    }
  }


void init_double_array(int size, double array[]){
    
    for(int i = 0; i < size; i++){
        array[i] = 0.0;
    }
}