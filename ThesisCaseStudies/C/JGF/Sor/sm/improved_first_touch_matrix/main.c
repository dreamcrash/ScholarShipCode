/* 
 * File:   main.c
 * Author: Dreamcrash
 *
 * Created on 20 of October of 2016, 15:57
 *
 * NOTE : The SOR parallelized versions (e.g., SM and DM) that are based on 
 * this sequential version, can only use a static block distribution by threads/processes. 
 * 
 * This distribution should distribute a unique contiguous block per thread/processes 
 * in a thread/process ID ascendency order. Thus, a distribution equivalent to the one 
 * of choosing the OpenMP 'parallel for static' with the a default chunk size. 
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "Sor.h"


int main(int argc, char** argv) {

    // If the user wants to validate the simulation 
    int validation  = (argc > 1) ? (strcmp(argv[1],"-v") == 0) :1;	
    int size        = (argc > 2) ? atoi(argv[2]) : 0;   // dim problem
    int nthreads    = (argc > 3) ? atoi(argv[3]) : 2;
    run(size, validation, nthreads);
    
    return (EXIT_SUCCESS);
}

