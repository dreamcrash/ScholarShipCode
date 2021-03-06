/** 
 *  File:   main.c
 *  Code Produced by Dreamcrash       
 *                                                                         
 *                   at                                       
 *                                                                         
 *  Uminho University. Using as base the code produced by    
 *  Java Grande Benchmarking Project
 * 
 * Molecular Dynamic (MD) using: 
 * -> Structure of Arrays (SoA) to represent the particles;
 * -> Optimizations:
 *      -> 3 Newton's law;
 * 	-> Use of Arrays to represent the particles (cache locality);
 * 
 * -> Distribution Memory :
 * 	-> Each Process will performance the entire MD simulation, with the
 * force computation iterations being split across processes.
 * 	-> After the end of the force computation: 
 * 		-> the variables {fx, fy, fz} of the Particles;
 * 		-> and the variables {epot, vir, interactions} of the MD Object;
 *              will be reduced across all Processes (All_reduce).
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <omp.h>
#include <mpi.h>
#include "Structs.h"
#include "main.h"
#include "MDSoA.h"
#include "functions_mpi.h"


/** Data sizes:
 *      8  -> 2048   particles;
 *  	13 -> 8788   particles;
 * 	17 -> 19652  particles;
 *  	40 -> 256000 particles;
 *  	50 -> 500000 particles;
 * Arguments 
 * 	argv[1] -> Validation (or not) of MD simulation
 * 	argv[2] -> Simulation Size
 * 	argv[3] -> # of Iterations
 **/

int main(int argc, char** argv){
    
    MPI_Init(&argc, &argv);
    
    int datasizes []    = {8,13,17,40,50};
    int validation  	= (argc > 1) ? (strcmp(argv[1],"-v") == 0) :1;
    int size 	  	= (argc > 2) ? atoi(argv[2]) : 0;
    int iterations 	= (argc > 3) ? atoi(argv[3]) : 50;
    
    MD *md = newMD(size, datasizes,iterations);
    
    /** Running simulation */
    MPI_Barrier(MPI_COMM_WORLD);
    double start = omp_get_wtime();
    runMD (md);
    
    MPI_Barrier(MPI_COMM_WORLD);
    double end = omp_get_wtime();
    
    if(md->processID == 0) 
    {
        printf("Time:%f\n",end-start);
        
        if(validation)
        {
            validate(md);
        }
    }
    
    freeMD(md);
  
    MPI_Finalize();

    return (EXIT_SUCCESS);
}

