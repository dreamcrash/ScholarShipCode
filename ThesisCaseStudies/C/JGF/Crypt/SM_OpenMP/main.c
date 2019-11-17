/* 
 * File:   main.c
 * Author: Bruno Medeiros
 *
 * Created on 14 de Outubro de 2016, 17:22
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
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

int getNumberOfThreadsToRunInParallel(int argc, char** argv)
{
    return (argc > 3) ? atoi(argv[3]) : 2;
}

int main(int argc, char** argv) {

    int validation          = shouldExecutionBeValidated(argc, argv);	
    int size                = getInputSize(argc, argv);
    int numThreads          = getNumberOfThreadsToRunInParallel(argc, argv);
    
    run(size, validation, numThreads);
    
    return (EXIT_SUCCESS);
}

