/* 
 * File:   SparseMatmult.h
 * Author: Dreamcrash
 *
 * Created on 16 de October de 2016, 23:40
 */

#ifndef SPARSEMATMULT_H
#define SPARSEMATMULT_H

#ifdef __cplusplus
extern "C" {
#endif

typedef struct SparseMatmult
{
    int datasizes_M;  
    int datasizes_N;
    int datasizes_nz;
    int SPARSE_NUM_ITER;
    long RANDOM_SEED;
    
    double *x;
    double *y;
    double *val;
    int *col;
    int *row;
    int size;
    
    // MPI CCC's
    double  *buf_val;
    double  *p_y;
    int     *buf_col;
    int     *buf_row;
    
    int p_datasizes_nz;
    int ref_p_datasizes_nz;
    int rem_p_datasizes_nz;
    double ytotal;
     
}Sparse;    
    
    
    
void    run                 (int size, int validation);
void    free_sparse         (Sparse **sparse, int rank);
int     JGFinitialise       (Sparse *sparse, int rank, int nprocess);
double *RandomVector        (int N);
int     abs_                (int value);
double  abs_double          (double value);
void    JGFkernel           (Sparse *sparse, int rank);
double  test                (double y[], double val[], int row[],
                                int col[], double x[], int NUM_ITERATIONS, 
                                int var_length, int buf_row[], double p_y[], 
                                int buf_row_length, int y_length, int rank);
void    JGFvalidate         (Sparse *sparse);



#ifdef __cplusplus
}
#endif

#endif /* SPARSEMATMULT_H */

