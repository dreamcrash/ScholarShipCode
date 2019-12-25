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
    double ytotal;
    
}Sparse;    
    
    
    
void    run             (int size, int validation);
void    free_sparse     (Sparse **sparse);
int     JGFinitialise   (Sparse *sparse);
double *RandomVector    (int N);
int     abs_            (int value);
double  abs_double      (double value);
void    JGFkernel       (Sparse *sparse);
double  test            (double y[], double val[], int row[],
                            int col[], double x[], int NUM_ITERATIONS, 
                            int var_length);
void    JGFvalidate         (Sparse *sparse);
void    init_double_array   (int size, double array[]);



#ifdef __cplusplus
}
#endif

#endif /* SPARSEMATMULT_H */

