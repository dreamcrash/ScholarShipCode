/** 
 *  File:   main.c
 *  Code Produced by dreamcrash       
 *                                                                         
 *                   at                                       
 *                                                                         
 *  Uminho University. Using as base the code produced by    
 *  Java Grande Benchmarking Project
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "mpi.h"
#include "main.h"
#include "SeriesTest.h"

int main(int argc, char** argv) {
    
    MPI_Init(&argc, &argv); // MPI ccc's
    
    // If the user wants to validate the simulation 
    int validation  	= (argc > 1) ? (strcmp(argv[1],"-v") == 0) :1;	
    int size 	  	= (argc > 2) ? atoi(argv[2]) : 0;   // dim problem
    const int datasizes[]={10000,100000,1000000, 2000000, 2500000};
    
    run(datasizes[size], validation);
    
    MPI_Finalize();
    return (EXIT_SUCCESS);
}

