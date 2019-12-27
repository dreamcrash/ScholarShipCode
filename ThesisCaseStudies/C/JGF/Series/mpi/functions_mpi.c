#include "mpi.h"
#include "functions_mpi.h"
#include <stdio.h>

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

// finding the master process
int master(){
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    return (rank == 0);
}

void reduceTestArray(const int process, const int size, double TestArray [2][size]){
    
    /**
     * When the communicator is an intracommunicator, you can perform 
     * a reduce operation in-place (the output buffer is used as the input buffer). 
     * Use the variable MPI_IN_PLACE as the value of the root process sendbuf. 
     * In this case, the input data is taken at the root from the receive buffer, 
     * where it will be replaced by the output data.
     */
    // Use the variable MPI_IN_PLACE as the value of the root process sendbuf.
    
    if(process == 0)
    {
        MPI_Reduce(MPI_IN_PLACE, TestArray, size * 2, 
                MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD); 
    }
    else
    {
        MPI_Reduce(TestArray, TestArray, size * 2, 
                MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD); 
    }
    

}