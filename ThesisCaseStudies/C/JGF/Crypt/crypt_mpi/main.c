/* 
 * File:   main.c
 * Author: Bruno Medeiros.
 *
 * Created on 14 de Outubro de 2016, 17:22
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>
#include "Crypt.h"

/**
 * 
 * @author Bruno Medeiros
 * 
 * This code is the result of reusing the Crypt case study from 
 * JGF JAVA Benchmarks and adapting it to C.
 */

int shouldExecutionBeValidated(int argc, char** argv)
{
    return (argc > 1) ? (strcmp(argv[1],"-v") == 0) :1;
}

int getInputSize(int argc, char** argv)
{
    return (argc > 2) ? atoi(argv[2]) : 0;
}

int main(int argc, char** argv) {

    /* Initialize MPI */
    MPI_Init(&argc, &argv);
    
    int validation          = shouldExecutionBeValidated(argc, argv);	
    int size                = getInputSize(argc, argv);
    
    run(size, validation);
    
        /* Finalize MPI */
     MPI_Finalize();
    
    return (EXIT_SUCCESS);
}

