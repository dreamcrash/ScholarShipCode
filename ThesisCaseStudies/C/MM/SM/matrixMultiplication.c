
/* 
 * File:   matrixMultiplication.c
 * Author: Dreamcrash
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include "matrixMultiplication.h"

/*
        How does the matrix multiplication work:
 * 
 * 
 *      A = [a11 a12 a13]       B = [b11 b12]
 *          [a21 a22 a23]         = [b21 b22]
 *          [a31 a32 a33]         = [b31 b32]
 *          [a41 a42 a43]         
 *  
 *      A (4 * 3)               B (3 * 2)
 * 
 *      C = 4 * 2
 *      C = # rows of A * # col of B
 * 
 *      C = [a11 * b11 + a12 * b21 + a13 * b31 + a11 * b12 + a12 * b22 + a13 * b32]
 *          [a21 * b11 + a22 * b21 + a23 * b31 + a21 * b12 + a22 * b22 + a23  *b32]
 *          [....]
 *          [....]
 * 
 *      You should note that in order to do A * B, the # of columns in A must be equal to the # of rows in B
 * 
 */
 
/**
 * Set matrix to zeros
 * @param rows    : The number of rows of the matrix
 * @param cols    : The number of columns of the matrix
 * @param matrix       : The pointer matrix
 */
void cleanMatrix(const int rows, const int cols, double matrix[rows][cols]){
    
    for(int i = 0; i < rows; i++)
        for(int j = 0; j < cols; j++)
            matrix[i][j] = 0;
}

/**
 * Copy the contend of one matrix to another
 * 
 * @param rows      : Number of lines to be copy
 * @param cols      : The number of elements per line
 * @param source    : The matrix from the values will be copied
 * @param dest      : The matrix where the values will be saved
 */
void copyMatrix(const int rows, const int cols, double source[rows][cols], double dest[rows][cols]){
    
    for(int i = 0; i < rows; i++)
        for(int j = 0; j < cols; j++)
            dest[i][j] = source[i][j];
}
	       
/**
 * Checking if two matrices are equal
 * @param rows  : Number of lines to be copy
 * @param cols  : The number of elements per line
 * @param a     : First matrix to be compare
 * @param b     : Second matrix to be compare
 * @return      : 0 if the matrices are not equal, 1 otherwise
 */
int equalsMatrix(  const int rows, const int cols, double a[rows][cols], double b[rows][cols]){
    
    for(int i = 0; i < rows; i++)
        for(int j = 0; j < cols; j++)
            if(a[i][j] != b[i][j])
                return 0;
    return 1;
}
			
/**
 * Sum all elements of a given matrix
 * @param rows  : The number of rows of the matrix
 * @param cols  : The number of columns of the matrix
 * @param m     : The pointer to the matrix
 * @return      : The result from adding all the elements of the matrix
 */
int sumMatrixElements(  const int rows, const int cols, double m[rows][cols]){
    
    int acc = 0;
    for(int i = 0; i < rows; i++)
        for(int j = 0; j < cols; j++)
            acc += m[i][j];
    return acc;
}

/** Matrix Multiplication algorithm **/
	    
/***********************  NAIVE VERSIONS  ***************************/

/**
 * Sequential matrix multiplication naive version (without any packing)
 * @param maxRowA   : The number of line of matrix A
 * @param maxColB   : The number of columns of matrix B
 * @param maxColA   : The number of columns of matrix A
 * @param A         : Matrix A
 * @param B         : Matrix B
 * @param C         : Matrix C
 */
void naiveMM(   const int maxRowA, const int maxColB, const int maxColA,
                double A[maxRowA][maxColA], double B[maxColA][maxColB], double C[maxRowA][maxColB]){
   
    for (int i = 0; i < maxRowA; i ++)
        for (int j = 0; j < maxColB; j++)
            for (int k = 0; k < maxColA; k++) 
                C[i][j] += A[i][k] * B[k][j];
}

void matrixMuliplication(   const int version,
                            const int maxRowA, const int maxColA, const int maxRowB, const int maxColB,
                            double A[maxRowA][maxColA], double B[maxColA][maxColB], double C[maxRowA][maxColB]){
    
    double begin = omp_get_wtime();
    if(version == 0)
    {
        naiveMM(maxRowA, maxColB, maxColA, A, B, C);
    }
    else if(version == 1)   
    {
        optimizedMM(maxRowB, maxColB, maxColA, A, B, C);
    }
    
    double end = omp_get_wtime();
    printf("Time: %f\n",end-begin);
}



/**********  OPTIMIZED Sequential VERSION **********/


/**
 *  Cache L3 packing, loads matrix B to a sub-matrix (bb2) in a continuously manner.
 * 
 * @param jj        : Columns begin
 * @param bb2       : the packing matrix where the matrix B will be loaded
 * @param maxRowB   : max number of rows of matrix B
 * @param maxColB   : max number of cols of matrix B
 * @param B         : the matrix from the values will be load
 */
void packingCacheL3(const  int jj, const  int maxRowB, const int maxColB, 
                     double bb2[maxRowB/4][tilej*4+1], double B[maxRowB][maxColB]){
            
    #pragma omp for schedule(dynamic)
    for(int k = 0; k < maxRowB; k += 4)
    {
        for(int j = 0; j < tilej; j++) 
        {
            bb2[k/4][j]         =   B[k]  [jj+j];
            bb2[k/4][j+tilej]	=   B[k+1][jj+j];
            bb2[k/4][j+tilej*2]	=   B[k+2][jj+j];
	    bb2[k/4][j+tilej*3]	=   B[k+3][jj+j];
        }
    }
}

void initBlock(double cc [tilei][tilej+1]){
    
    for(int i = 0; i<tilei; i++)      // Initialization of a tile 32 * 256 
        for(int j = 0; j<tilej; j++)  // that fits in L2 cache
            cc[i][j]=0;
}

/**
 *  Copy the results from the tile to the matrix C
 * @param ii        : the begin line from the matrix C
 * @param jj        : the begin column from the matrix C
 * @param cc        : the matrix with partial results
 * @param maxColB   : the max number of cols of matrix C (the as same as matrix B)
 */
void updateCCmatrix(const  int ii, const  int jj, const int maxColB, 
                    double cc [tilei][tilej+1],double C[][maxColB]){
    
    for(int i = 0; i < tilei; i++)
        for(int j = 0; j < tilej; j++)
            C[ii+i][j+jj] += cc[i][j];       
}
			
/**
 *  This micro Tiling will perform  a tile (4x2) over the original tile (256 x 32)
 * @param ii 	: the begin line from the matrix C
 * @param bb2	: the matrix packing from matrix B
 * @param cc	: the matrix packing from matrix C
 */
void microTiling(const int ii, const int maxRowB, const int maxRowA, const int maxColA, 
                 double bb2 [maxRowB/4][tilej*4+1], double cc[tilei][tilej+1], 
                 double A[maxRowA][maxColA]){
    
    for (int k = 0; k < maxRowB; k += 4)            // Loads 4 Columns			
    {
	double *bk = bb2[k/4];                  
        for (int i = ii; i < ii + tilei; i += 2)    // and 2 lines
        {
             // Loading a tile 2 x 4 to the registers
            double a11  =   A[i][k];					
	    double a12  =   A[i][k+1];
	    double a13  =   A[i][k+2];
	    double a14  =   A[i][k+3];
	               
            double a21  =   A[i+1][k];
            double a22  =   A[i+1][k+1];
            double a23  =   A[i+1][k+2];
            double a24  =   A[i+1][k+3];
	                    
	    double *ci  = cc[i-ii];
	    double *ci1 = cc[i-ii+1];
	                       
	    // performing the multiplication of the micro-tile
            for (int j = 0; j <  tilej + 1; j++)	
            { 
                // Load 4 lines of tile B to L1 cache
                double bkj  = bk[j];
                double bk1j = bk[j+1 * tilej];
                double bk2j = bk[j+2 * tilej];
                double bk3j = bk[j+3 * tilej];
	                    
                double aux1 = ci[j];
                aux1 += a11 * bkj;
                aux1 += a12 * bk1j;
                aux1 += a13 * bk2j;
                aux1 += a14 * bk3j;
                ci[j] = aux1;
                
	        double aux2 = ci1[j];
                aux2 += a21 * bkj;
                aux2 += a22 * bk1j;
                aux2 += a23 * bk2j;
                aux2 += a24 * bk3j;
                ci1[j] = aux2;
            } // for j	                   
        } // for i
    } // for k
}
		
void optimizedMM( const int maxRowA, const int maxColB, const int maxColA, 
                  double A[maxRowA][maxColA], double B[maxColA][maxColB], double C[maxRowA][maxColB]){

    const int maxColC = maxColB, maxRowC = maxRowA;
    const int maxRowB = maxColA;
    
    // Packing matrices (i.e., cc && bb) to hold a chunk of the values of matrix C and B.
    double (*bb) [tilej*4+1] = malloc(sizeof *bb * (maxRowB/4)); 
   
    // Permutes horizontally the tile of matrix C
    #pragma omp parallel
    {   
        double (*cc) [tilej+1] = malloc(sizeof *cc * tilei); 
        for(int jj = 0; jj < maxColC; jj += tilej)     
        {
            // Cache L3 packing, loads matrix B to bb2 in a continuously manner
            packingCacheL3(jj, maxRowB, maxColB, bb, B);    

            #pragma omp for schedule(dynamic)
            for(int ii = 0; ii < maxRowC; ii += tilei) 
            {
                initBlock       (cc);
                microTiling     (ii, maxRowB, maxRowA, maxColA, bb, cc, A);  
                updateCCmatrix  (ii, jj, maxColB, cc, C);               
            }
        }
        free(cc);
    }
    free(bb);
}