/* 
 * File:   main.c
 * Author: Dreamcrash
 *  
 * Matrix multiplication versions with two versions :
 * => Version 0 -> naive non-optimized version;
 * => Version 1 -> Shared Memory + Distributed Memory version with tiling version:
 *                  -> Matrices A and C are divided in chunks of lines among processes;
 *                  -> Every process has the entire matrix B.
 */

#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <omp.h>
#include "main.h"
#include "processMPI.h"
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

    MPI_Init(&argc, &argv);
    
    // User options
    const int validate      = (argc > 1) ? atoi(argv[1]) : 0;
    const int version       = (argc > 2) ? atoi(argv[2]) : 1;
    int maxRowA             = (argc > 3) ? atoi(argv[3]) : 1024;
    int maxColA             = (argc > 4) ? atoi(argv[4]) : 1024;
    int maxRowB             = (argc > 5) ? atoi(argv[5]) : 1024;
    int maxColB             = (argc > 6) ? atoi(argv[6]) : 1024;
    const int numThreads = (argc > 7) ? atoi(argv[7]) : 2;
    int maxRowC             = maxRowA;
    int maxColC             = maxColB;
    const int idProcess     = getProcessId();
    const int numProcess    = numberProcess();
    
    /** 
     * Regardless of the version the master process will always have the entire
     * matrices. However, in the DM version the slaves will only have a sub-set of
     * the matrix A and C. Those matrices will be divided by lines.
     */
    if(version == 1 && idProcess != 0)
    {
        // Calculate the subset size of the matrices to be split
        maxRowC = maxRowA = myMatrixSize(maxRowA/tilei, idProcess, numProcess) * tilei;
        // If this slave will not compute any row of the matrix 
        // (happens when there is more processes than chunk of worker to be performed)
        if(maxRowC == 0) maxColA = maxColC = 0;
    }
    
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
        if(idProcess == 0)
        {
            fillupRandomly (MIN, MAX, RANDOM_SEED, maxRowA, maxColA, A);
            fillupRandomly (MIN, MAX, RANDOM_SEED, maxRowB, maxColB, B);
        }
        else
        {
            cleanMatrix   (maxRowA, maxColA, A); 
            cleanMatrix   (maxRowB, maxColB, B);
        }
        
        cleanMatrix   (maxRowA, maxColB, C);
        omp_set_dynamic(0);              /** Explicitly disable dynamic teams **/
        omp_set_num_threads(numThreads); /** Use N threads for all parallel regions **/
        matrixMuliplication(version, idProcess, numProcess, maxRowA, maxColA, maxRowB, maxColB, A, B, C);
        
        if(validate && idProcess == 0)
        {
            validateMM(maxRowA, maxColB, maxColA, A, B, C);    
        }
        
        free(A);
        free(B);
        free(C);
    }
    
    MPI_Finalize();
     
    return (EXIT_SUCCESS);
}

