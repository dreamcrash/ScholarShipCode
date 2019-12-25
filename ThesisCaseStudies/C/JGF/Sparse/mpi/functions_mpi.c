/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

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


void send_data_to_slaves(int *row, int *col, double *val, int dim ,int slaveID)
{    
    MPI_Send(row, dim, MPI_INT,     slaveID, slaveID, MPI_COMM_WORLD);
    MPI_Send(col, dim, MPI_INT,     slaveID, slaveID, MPI_COMM_WORLD);
    MPI_Send(val, dim, MPI_DOUBLE,  slaveID, slaveID, MPI_COMM_WORLD);
}

void recv_data_from_master(int *row, int *col, double *val, int dim , 
        int slaveID)
{
    
    MPI_Recv(row, dim, MPI_INT, 0, 
            slaveID, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    
    MPI_Recv(col, dim, MPI_INT, 0, 
            slaveID, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    
    MPI_Recv(val, dim, MPI_DOUBLE, 0, 
            slaveID, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
}