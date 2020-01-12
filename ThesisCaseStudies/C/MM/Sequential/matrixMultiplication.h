
/* 
 * File:   matrixMultiplication.h
 * Author: Dreamcrash
 *
 */

#ifndef MATRIX_H
#define MATRIX_H

#ifdef __cplusplus
extern "C" {
#endif
    
#define tilei 32
#define tilej 256    
    
void cleanMatrix        (const int rows, const int cols, double matrix[rows][cols]);
void copyMatrix         (const int rows, const int cols, double source[rows][cols], double dest[rows][cols]);
int equalsMatrix        (const int rows, const int cols, double a[rows][cols], double b[rows][cols]);
int sumMatrixElements   (const int rows, const int cols, double m[rows][cols]);
void naiveMM            (   const int maxRowA, const int maxColB, const int maxColA,
                            double A[maxRowA][maxColA], double B[maxColA][maxColB], double C[maxRowA][maxColB]);
void matrixMuliplication(const int version, const int maxRowA, const int maxColA, const int maxRowB, const int maxColB,
                            double A[maxRowA][maxColA], double B[maxRowB][maxColB], double C[maxRowA][maxColB]);
void packingCacheL3     (const int jj, const int maxRowB, const int maxColB, 
                            double bb2[maxRowB/4][tilej*4+1], double B[maxRowB][maxColB]);
void initBlock          (double cc [tilei][tilej+1]);
void updateCCmatrix     (const int ii, const int jj, const int maxColB, double cc [tilej][tilej+1], double C[][maxColB]);
void microTiling        (const int ii, const int maxRowB, const int maxRowA, const int maxColA, 
                            double bb2 [maxRowB/4][tilej*4+1], double cc[tilei][tilej+1], double A[maxRowA][maxColA]);
void optimizedMM        ( const int maxRowA, const int maxColB, const int maxColA, 
                            double A[maxRowA][maxColA], double B[maxColA][maxColB], double C[maxRowA][maxColB]);

#ifdef __cplusplus
}
#endif

#endif /* MATRIX_H */

