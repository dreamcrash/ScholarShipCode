/* 
 * File:   main.c
 * Author: Dreamcrash
 *  
 * Matrix multiplication versions with two versions :
 * => Version 0 -> naive non-optimized version;
 * => Version 1 -> Shared Memory with tiling version.
 */

#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include "main.h"
#include "matrixMultiplication.h"
#include "RandomMatrix.h"
#include "main.h"

void validateMM(const int maxRowA, const int maxColB, const int maxColA,
                double A[maxRowA][maxColA], double B[maxColA][maxColB], double C[maxRowA][maxColB]) {

    double (*validationMatrixC) [maxColB] = malloc(sizeof *validationMatrixC * maxRowA);

    cleanMatrix(maxRowA, maxColB, validationMatrixC);
    naiveMM(maxRowA, maxColB, maxColA, A, B, validationMatrixC);

    if (equalsMatrix(maxRowA, maxColB, C, validationMatrixC)) 
    {
        printf("Valid Matrix Multiplication => {%d}\n", sumMatrixElements(maxRowA, maxColB, C));
    } 
    else 
    {
        printf("Invalid\n");
    }
}

int main(int argc, char** argv) {

    // User options
    const int validate      = (argc > 1) ? atoi(argv[1]) : 0;
    const int version       = (argc > 2) ? atoi(argv[2]) : 1;
    const int maxRowA       = (argc > 3) ? atoi(argv[3]) : 1024;
    const int maxColA       = (argc > 4) ? atoi(argv[4]) : 1024;
    const int maxRowB       = (argc > 5) ? atoi(argv[5]) : 1024;
    const int maxColB       = (argc > 6) ? atoi(argv[6]) : 1024;
    const int numThreads = (argc > 7) ? atoi(argv[7]) : 2;
    
    // Matrices allocation
    double (*A) [maxColA] = malloc(sizeof *A * maxRowA);
    double (*B) [maxColB] = malloc(sizeof *B * maxRowB);
    double (*C) [maxColB] = malloc(sizeof *C * maxRowA);
    
    // Dealing with fail memory allocation
    if(!A || !B || !C)
    {
        if(A)   free(A);
        if(B)   free(B);
        if(C)   free(C);
        return (EXIT_FAILURE);
    }
    else
    {
        fillupRandomly (MIN, MAX, RANDOM_SEED, maxRowA, maxColA, A);
        fillupRandomly (MIN, MAX, RANDOM_SEED, maxRowB, maxColB, B);
        cleanMatrix   (maxRowA, maxColB, C);
       
        omp_set_dynamic(0);              /** Explicitly disable dynamic teams **/
        omp_set_num_threads(numThreads); /** Use N threads for all parallel regions **/
        matrixMuliplication(version, numThreads, maxRowA, maxColA, maxRowB, maxColB, A, B, C);
        
        if(validate)
        {
            validateMM(maxRowA, maxColB, maxColA, A, B, C);    
        }
        
        free(A);
        free(B);
        free(C);
    }
    
    return (EXIT_SUCCESS);
}

