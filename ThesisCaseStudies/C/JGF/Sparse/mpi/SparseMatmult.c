#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <mpi.h>
#include <time.h>
#include "functions_mpi.h"
#include "SparseMatmult.h"




void run (int size, int validation){
        
    int rank        = getProcessId();
    int nprocess    = numberProcess();
    Sparse *sparse = malloc(sizeof(Sparse));
    sparse->size = size;
    
    if(!JGFinitialise(sparse, rank, nprocess))
    {
        printf("Memory problems at JGFinitialise! \n");
    }
    else
    {
        JGFkernel(sparse, rank);

        if(rank == 0 && validation)
        {
            JGFvalidate(sparse);
        }
    
        free_sparse(&sparse, rank);
    }
   
}

void free_sparse(Sparse **sparse, int rank){
    
    if(rank == 0)
    {
        free((*sparse)->buf_col);
        free((*sparse)->buf_row);
        free((*sparse)->buf_val);
    }
    free((*sparse)->col);
    free((*sparse)->row);
    free((*sparse)->val);
    free((*sparse)->x);
    free((*sparse)->y);
    free((*sparse)->p_y);
    free((*sparse));
    *sparse         = NULL;
}

int JGFinitialise(Sparse *sparse, int rank, int nprocess){

    int datasizes_M[]   = { 50000, 100000,  500000, 1000000, 1500000 };  
    int datasizes_N[]   = { 50000, 100000,  500000, 1000000, 1500000 };
    int datasizes_nz[]  = {250000, 500000, 2500000, 5000000, 7500000 };
   
    int size                = sparse->size;
    
    sparse->datasizes_M     = datasizes_M[size];
    sparse->datasizes_N     = datasizes_N[size];
    sparse->datasizes_nz    = datasizes_nz[size];
    
    // MPI CCC's
    sparse->p_datasizes_nz      = (sparse->datasizes_nz + nprocess -1) /nprocess;
    sparse->ref_p_datasizes_nz  = sparse->p_datasizes_nz;
    sparse->rem_p_datasizes_nz  = sparse->p_datasizes_nz - 
            ((sparse->p_datasizes_nz * nprocess) - sparse->datasizes_nz);
    
    if(rank==(nprocess-1) 
            && (sparse->p_datasizes_nz *(rank+1) > sparse->datasizes_nz))
    {
         sparse->p_datasizes_nz = sparse->rem_p_datasizes_nz;
    }
    
    sparse->SPARSE_NUM_ITER = 200;
    sparse->RANDOM_SEED     = 10101010;
    
    srand(sparse->RANDOM_SEED);
   
    sparse->x       = RandomVector(sparse->datasizes_N);
    sparse->y       = malloc(sizeof(double) * sparse->datasizes_M);
    sparse->p_y     = malloc(sizeof(double) * sparse->datasizes_M);
    
    for(int i = 0; i < sparse->datasizes_M; i++)
    {
        sparse->y[i] = sparse->p_y[i] = 0.0;
    }
    
    sparse->val     = malloc(sizeof(double) * sparse->p_datasizes_nz);
    sparse->col     = malloc(sizeof(int)    * sparse->p_datasizes_nz);
    sparse->row     = malloc(sizeof(int)    * sparse->p_datasizes_nz);
    
    for(int i = 0; i < sparse->p_datasizes_nz; i++)
    {
        sparse->val[i] = sparse->col[i] = sparse->row[i] = 0;
    }
    
    if(!(sparse->x  && sparse->y && sparse->p_y
            && sparse->val && sparse->col && sparse->row))
    {
        return 0;
    }
    
    
    
    if(rank==0) 
    {
        sparse->buf_val = malloc(sizeof(double) * sparse->datasizes_nz);
        sparse->buf_col = malloc(sizeof(int)    * sparse->datasizes_nz);
        sparse->buf_row = malloc(sizeof(int)    * sparse->datasizes_nz);
        
        if(!(sparse->buf_val && sparse->buf_col && sparse->buf_row)) return 0;
        
        for(int i = 0; i < sparse->datasizes_nz; i++)
        {
            sparse->buf_val[i] = sparse->buf_col[i] = sparse->buf_row[i] = 0;
        }
        
        int p_datasizes_nz = sparse->p_datasizes_nz;
        for(int slaveID = 0; slaveID < nprocess; slaveID++) 
        {
            if(slaveID == nprocess-1) 
            {
                p_datasizes_nz = sparse->rem_p_datasizes_nz;
            } 
            
            for (int i = 0; i< p_datasizes_nz; i++) 
            {
                sparse->buf_row[i + slaveID * sparse->ref_p_datasizes_nz] = 
                        sparse->row[i] = abs_(rand()) % sparse->datasizes_M;
                
                sparse->buf_col[i + slaveID * sparse->ref_p_datasizes_nz] = 
                        sparse->col[i] = abs_(rand()) % sparse->datasizes_N;
                
                sparse->buf_val[i + slaveID * sparse->ref_p_datasizes_nz] = 
                        sparse->val[i] = (((double)rand()/(double)RAND_MAX));
         
            }
            // slaves
            if(slaveID != 0)
            {
                send_data_to_slaves(sparse->row, sparse->col, sparse->val, 
                        p_datasizes_nz, slaveID);
            }
        }
        p_datasizes_nz = sparse->ref_p_datasizes_nz;
        
        // master generation
        for (int i = 0; i < p_datasizes_nz; i++) 
        {
            sparse->row[i] = sparse->buf_row[i]; 
            sparse->col[i] = sparse->buf_col[i];
            sparse->val[i] = sparse->buf_val[i];
        }	
    }
    else
    {
        sparse->buf_val = sparse->buf_col = sparse->buf_row = NULL;
        recv_data_from_master  (sparse->row, sparse->col, sparse->val, 
                                sparse->p_datasizes_nz, rank);
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

void JGFkernel(Sparse *sparse, int rank){
    
    sparse->ytotal = test(sparse->y, sparse->val, sparse->row,
                         sparse->col, sparse->x, sparse->SPARSE_NUM_ITER,
                         sparse->p_datasizes_nz, sparse->buf_row, sparse->p_y,
                         sparse->datasizes_nz, sparse->datasizes_M, rank);     
}

double test (   double y[], double val[], int row[],
                int col[], double x[], int NUM_ITERATIONS, 
                int var_length, int buf_row [], double p_y[], 
                int buf_row_length, int y_length, int rank){
    
    int nz = var_length;
    double ytotal = 0;
    
    MPI_Barrier(MPI_COMM_WORLD);
    const double start  = omp_get_wtime();
    
    for (int reps = 0; reps < NUM_ITERATIONS; reps++)
    {
        for (int i = 0; i < nz; i++)
        {
            p_y[ row[i] ] += x[ col[i] ] * val[i];
        }
        MPI_Allreduce(p_y, y, y_length, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    }
    
    
    MPI_Barrier(MPI_COMM_WORLD);
    const double end    = omp_get_wtime();
   

    if(rank == 0)
    { 
        printf("Timer:%f\n", (end - start));
        for (int i = 0; i < buf_row_length; i++) 
        {
            ytotal += y[  buf_row[i] ];
        }
    }
    return ytotal;
}

void JGFvalidate(Sparse *sparse){
    
    double refval[] = { 75.205209621081394, 150.325289707762209, 
                        750.707296942218136, 1499.080290091720144, 
                        2250.518217459347852};
    
    double dev = abs_double(sparse->ytotal - refval[sparse->size]);
    if (dev > 1.0e-9)
    {
        printf("Validation failed \n");
	printf("ytotal = %.15f %.15f %d \n", sparse->ytotal, dev, sparse->size);
    }
  }