/** 
 *  File:   main.c
 *  Code Produced by Dreamcrash.        
 *                                                                         
 *                   at                                       
 *                                                                         
 *  Uminho University. Using as base the code produced by    
 *  Java Grande Benchmarking Project
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "main.h"
#include "Linpack.h"

int main(int argc, char** argv) {

    // If the user wants to validate the simulation 
    int validation          = (argc > 1) ? (strcmp(argv[1],"-v") == 0) :1;	
    int size                = (argc > 2) ? atoi(argv[2]) : 0;   // dim problem
    
    run(size, validation);

    return (EXIT_SUCCESS);
}

