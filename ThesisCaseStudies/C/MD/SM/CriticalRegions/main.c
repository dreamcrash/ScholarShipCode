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
 *      -> critical regions to deal with the race conditions;
 *      -> Dividing the force calculation between particles among threads in a static fashion.
 * -> Optimizations:
 *      -> 3 Newton's law;
 * 	-> Use of Arrays to represent the particles (cache locality);
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <omp.h>
#include "Structs.h"
#include "main.h"
#include "MDSoA.h"

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
    
    int datasizes []    = {8,13,17,40,50};
    int validation  	= (argc > 1) ? (strcmp(argv[1],"-v") == 0) :1;
    int size 	  	= (argc > 2) ? atoi(argv[2]) : 0;
    int iterations 	= (argc > 3) ? atoi(argv[3]) : 50;
    int numThreads      = (argc > 4) ? atoi(argv[4]) : 2;
    
    MD *md = newMD(size, datasizes,iterations, numThreads);
    
    /** Running simulation */
    double start = omp_get_wtime();
    runMD (md);
    double end = omp_get_wtime();
    
    printf("%f\n",end-start);
    
    if(validation)
    {
        validate(md);
    }
    
    freeMD(md);
    return (EXIT_SUCCESS);
}

