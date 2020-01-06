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
 * -> SM;
 *      -> Using data redundancy to deal with the race conditions;
 *          -> In this version the forces are reduced manually: 
 *              -> avoiding therefore the tuning of the thread stack size needed 
 *                  on the version with OpenMP 4.0+ features;
 *              -> Uses less memory since the master thread uses the original particles forces (i.e,  does not replicate them);
 *              -> The downside is that this version is more intrusive and less readable than the one using the OpenMP 4.0+ features. 
 *      -> Improved load balancing by dividing the force calculation between particles among threads in a dynamic fashion.
 * 
 * -> Distribution Memory :
 * 	-> Each Process will performance the entire MD simulation, with the
 * force computation iterations being split across processes.
 * 	-> After the end of the force computation: 
 * 		-> the variables {fx, fy, fz} of the Particles;
 * 		-> and the variables {epot, vir, interactions} of the MD Object;
 *              will be reduced across all Processes (All_reduce).
 * -> Optimizations:
 *      -> 3 Newton's law;
 * 	-> Use of Arrays to represent the particles (cache locality);
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
    int numThreads      = (argc > 4) ? atoi(argv[4]) : 2;
    
    MD *md = newMD(size, datasizes,iterations, numThreads);
    
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

