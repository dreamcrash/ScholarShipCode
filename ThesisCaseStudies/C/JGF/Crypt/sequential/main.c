/* 
 * File:   main.c
 * Author: dreamcrash
 *
 * Created on 14 de Outubro de 2016, 17:22
 */

/**
 * 
 * @author dreamcrash
 * 
 * This code is the result of reusing the Crypt case study from 
 * JGF JAVA Benchmarks and adapting it to C.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "Crypt.h"

int shouldExecutionBeValidated(int argc, char** argv)
{
    return (argc > 1) ? (strcmp(argv[1],"-v") == 0) :1;
}

int getInputSize(int argc, char** argv)
{
    return (argc > 2) ? atoi(argv[2]) : 0;
}

int main(int argc, char** argv) {

    int validation          = shouldExecutionBeValidated(argc, argv);	
    int size                = getInputSize(argc, argv);
    
    run(size, validation);
    
    return (EXIT_SUCCESS);
}
