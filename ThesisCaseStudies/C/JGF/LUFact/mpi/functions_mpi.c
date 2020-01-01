#include "mpi.h"
#include "functions_mpi.h"
#include <stdio.h>

// Get process id
int getProcessId(){
    int rank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    return rank;
}

// Get number of process
int numberProcess(){
    int numProc = 1;
    MPI_Comm_size(MPI_COMM_WORLD, &numProc);
    return numProc;
}

void sendAndRecvSubMatrixCycleManner (int max_row, int max_col, double **a){
    
    int numberProcesses = numberProcess();
    int idProcess       = getProcessId();
 
    if(idProcess == 0)
    {
        for(int slaveID = 1; slaveID < numberProcesses; slaveID++)		
        {	   
            for(int line = slaveID ; line < max_row; line += numberProcesses) 
            {	
                MPI_Send(a[line], max_col, MPI_DOUBLE, slaveID, slaveID, MPI_COMM_WORLD);	
            }	
        }
    }
    else
    {
        for(int line = idProcess; line < max_row; line += numberProcesses)
        {   
            MPI_Recv(a[line], max_col, MPI_DOUBLE, 0, idProcess, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
    }
}

void send_buffer_to_worker(double *buf_col_k, int n, int worker){
    MPI_Bcast(buf_col_k, n, MPI_DOUBLE, worker, MPI_COMM_WORLD);
}

void send_l_to_worker(int *tmp_l, int worker){
    MPI_Bcast(tmp_l, 1, MPI_INT, worker, MPI_COMM_WORLD);
}

void master_collecting_results(int totalProcess, int max_row, int max_col, double **a){

    for(int slaveID = 1; slaveID < totalProcess; slaveID++)
    {
        for(int line = slaveID; line < max_row; line += totalProcess)
        {
            MPI_Recv(a[line], max_col, MPI_DOUBLE, slaveID, slaveID, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	}
    } 
}

void slaves_send_results( int slaveID, int totalProcess, int max_row, int max_col, double **a){
    
    for(int line = slaveID ; line < max_row; line += totalProcess) 
    {	
        MPI_Send(a[line], max_col, MPI_DOUBLE, 0, slaveID, MPI_COMM_WORLD);	
    }
}