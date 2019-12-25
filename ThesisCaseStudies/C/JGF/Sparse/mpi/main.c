/* 
 * File:   main.c
 * Author: Dreamcrash
 *
 * Created on 16 de October de 2016, 23:38
 */

/**
 * 
 * @author Dreamcrash
 * 
 * This code is the result of reusing the Sparse case study from 
 * JGF JAVA Benchmarks and adapting it to C.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>
#include "SparseMatmult.h"

int main(int argc, char** argv) {

    /* Initialize MPI */
    MPI_Init(&argc, &argv);
    
    // If the user wants to validate the simulation 
    int validation  = (argc > 1) ? (strcmp(argv[1],"-v") == 0) :1;	
    int size        = (argc > 2) ? atoi(argv[2]) : 0;   // dim problem
    
    run(size, validation);
    
     /* Finalize MPI */
     MPI_Finalize();
    
    return (EXIT_SUCCESS);
}

