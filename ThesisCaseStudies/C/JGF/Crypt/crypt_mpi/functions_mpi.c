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

void scattering (Ideatest *ideatest, int rank, int numberProcesses){
    
    byte *p_plain1      = ideatest->p_plain1;
    int  p_array_rows   = ideatest->p_array_rows;
    
    // master process
    if(rank == 0)
    {
        byte *plain1    = ideatest->plain1;
        
        // Internal Master copy
        for(int i = 0; i < p_array_rows; i++)	
        {
            p_plain1[i] = plain1[i];
        }

        for(int slaveID = 1; slaveID < numberProcesses-1; slaveID++)
        {	   
            MPI_Send(   &plain1[slaveID * p_array_rows], p_array_rows, 
                        MPI_UNSIGNED_CHAR, 
                        slaveID, slaveID, MPI_COMM_WORLD);
        }
        
        if(numberProcesses > 1) // Dealing with the last process
        {
            MPI_Send(   &plain1[(numberProcesses-1) * p_array_rows], 
                        ideatest->rem_p_array_rows, 
                        MPI_UNSIGNED_CHAR, 
                        (numberProcesses-1), 
                        (numberProcesses-1), 
                        MPI_COMM_WORLD);
        }
    }
    else // Slaves
    {
        MPI_Recv(p_plain1, p_array_rows, MPI_UNSIGNED_CHAR, 
                    0, rank, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
 }

void gathering  (Ideatest *ideatest, int rank, int numberProcesses){
    
    byte *p_plain2      = ideatest->p_plain2;
    int p_array_rows    = ideatest->p_array_rows;
    
    // Master process
    if(rank == 0)
    {
        byte *plain2        = ideatest-> plain2;
        int p_plain2_length = ideatest->p_array_rows;
        
        for(int i = 0; i < p_plain2_length; i++) // Internal Master copy
	{
            plain2[i] = p_plain2[i];
	}
		
        for(int slaveID = 1; slaveID < numberProcesses-1; slaveID++)
        {	   
            MPI_Recv(&plain2[slaveID * p_array_rows], p_array_rows, 
                    MPI_UNSIGNED_CHAR, 
                    slaveID, slaveID, 
                    MPI_COMM_WORLD, 
                    MPI_STATUS_IGNORE);
        }
        
        if(numberProcesses > 1)
        {
            MPI_Recv(   &plain2[(numberProcesses-1) * p_array_rows], 
                        ideatest->rem_p_array_rows, 
                        MPI_UNSIGNED_CHAR, 
                        (numberProcesses-1), (numberProcesses-1), 
                        MPI_COMM_WORLD, 
                        MPI_STATUS_IGNORE);
	}
    }
    else // Slaves
    {
        MPI_Send (p_plain2, p_array_rows, MPI_UNSIGNED_CHAR, 
                    0, rank, MPI_COMM_WORLD);
    }
}
