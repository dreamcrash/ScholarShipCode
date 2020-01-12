#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include "processMPI.h"

int getProcessId(){
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    return rank;
}

int numberProcess(){
    int numProc;
    MPI_Comm_size(MPI_COMM_WORLD, &numProc);
    return numProc;
}

/**
 * Function used to calculate the size of the subset of a given matrix
 * @param totalSize         : The total number of lines in the matrix
 * @param processID         : The id of the current process
 * @param numProcesses      : The number of process running in parallel
 * @return                  : The number of lines that the sub-matrix should
 *                              have
 */
int myMatrixSize (const int totalSize, const int processID, const int numProcesses){
    
    if(totalSize >= numProcesses)
    {
        return	processID < totalSize % numProcesses ? totalSize / numProcesses + 1 : totalSize / numProcesses;
    }
    
    return (processID < totalSize) ? 1 : 0;
}

/**
 * Function which the master will split a given matrix across the slaves 
 * processes in a cycle manner by lines
 * 
 * This function takes advantage of the fact that the matrix in continuously 
 * allocated in memory, sending blocks of lines (tilling block) 
 * in one single msg instead of multiple msg per tilling block.
 * 
 * 
 * @param maxRow            : The number of rows of the matrix
 * @param maxCol            : The number of columns of the matrix
 * @param m                 : The pointer to the matrix
 * @param tile              : The tiling size to be used in the partition 
 * @param idProcess         : The id of the current process
 * @param numberProcesses   : The number of process running in parallel
 */
void masterSendSubMatricesToSlaves (const int maxRow, const int maxCol, double m[][maxCol], const int tile, const int processID, const int numProcesses){
    
    const int taskInc = numProcesses * tile; 
    const int msgSize = maxCol*tile; 
 
    if(processID == 0)
    {
        for(int slave = 1; slave < numProcesses; slave++)
            for(int j = slave * tile; j < maxRow; j+= taskInc)
                MPI_Send(m[j], msgSize, MPI_DOUBLE, slave, slave, MPI_COMM_WORLD); 
    }
    else
    {
        for(int i = 0; i < maxRow; i+= tile)
            MPI_Recv(m[i], msgSize, MPI_DOUBLE, 0, processID, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
}

/**
 * Function which the master will collect a given matrix across the slaves 
 * processes in a block manner
 * 
 * This function takes advantage of the fact that the matrix in continuously 
 * allocated in memory, sending blocks of lines (tilling block) 
 * in one single msg instead of multiple msg per tilling block.
 * 
 * @param maxRow            : The number of rows of the matrix
 * @param maxCol            : The number of columns of the matrix
 * @param m                 : The pointer to the matrix
 * @param tile              : The tiling size to be used in the partition 
 * @param processID         : The id of the current process
 * @param numProcesses   : The number of process running in parallel
 */
void masterCollectsSubMatrix(const int maxRow, const int maxCol, double m[][maxCol], const int tile, const int processID, const int numProcesses){
    
    const int taskInc = numProcesses * tile; 
    const int msgSize = maxCol*tile; 
    if(processID == 0)
    {
        for(int slave = 1; slave < numProcesses; slave++)
            for(int j = slave * tile; j < maxRow; j+= taskInc)
                MPI_Recv(m[j], msgSize, MPI_DOUBLE, slave, slave, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
    else
    {
        for(int j = 0; j < maxRow; j+= tile)
            MPI_Send(m[j], msgSize, MPI_DOUBLE, 0, processID, MPI_COMM_WORLD);
    }
}