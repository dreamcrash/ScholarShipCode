#include <stdio.h>
#include <mpi.h>
#include "functions_mpi.h"
#include "Crypt.h"


// Get process id
int getProcessId(){
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    return rank;
}

// Get number of process
int numberProcess(){
    int numProc;
    MPI_Comm_size(MPI_COMM_WORLD, &numProc);
    return numProc;
}

void scattering (Ideatest *ideatest, int numberProcesses){
    
    int sendcounts [numberProcesses];
    int	displs [numberProcesses];
    fillMPIsendCountsAndDispls( ideatest, sendcounts, displs, numberProcesses);
	
    // https://www.mpich.org/static/docs/v3.2.x/www3/MPI_Scatterv.html
    MPI_Scatterv (  ideatest->plain1, 
                    sendcounts, 
                    displs, 
                    MPI_UNSIGNED_CHAR,
                    ideatest->p_plain1,
                    ideatest->p_array_rows, 
                    MPI_UNSIGNED_CHAR, 
                    0,
                    MPI_COMM_WORLD
                );
 }

void gathering  (Ideatest *ideatest, int numberProcesses){
     
    int sendcounts [numberProcesses];
    int	displs [numberProcesses];
    
    fillMPIsendCountsAndDispls( ideatest, sendcounts, displs, numberProcesses);
    
    // https://www.mpich.org/static/docs/v3.2.x/www3/MPI_Gatherv.html
    
    MPI_Gatherv(    ideatest->p_plain2, 
                    ideatest->p_array_rows, 
                    MPI_UNSIGNED_CHAR, 
                    ideatest->plain2, 
                    sendcounts, 
                    displs, 
                    MPI_UNSIGNED_CHAR, 
                    0,
                    MPI_COMM_WORLD
                );
}

void fillMPIsendCountsAndDispls( Ideatest *ideatest, int sendcounts[], 
                                 int displs [], int numberProcesses){
    const int p_array_rows = ideatest->p_array_rows;
    for (int slaveID = 0; slaveID < numberProcesses; slaveID++) 
    {
	    sendcounts[slaveID] = p_array_rows;
	    displs[slaveID] = slaveID * p_array_rows;
    }
    // Likely the last process has a 
    // smaller chunk than the remaining processes
    sendcounts[numberProcesses-1] = ideatest->rem_p_array_rows;
}