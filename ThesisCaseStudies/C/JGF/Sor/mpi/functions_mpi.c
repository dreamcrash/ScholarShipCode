#include "mpi.h"
#include "functions_mpi.h"
#include "Sor.h"
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

void sendMatrixFromMasterToSlaves(  Sor* sor, const int rank, const int nprocess, 
                                    int p_G_lenght, int N, double p_G[p_G_lenght][N], 
                                    int M, double G[M][N]){
    if(rank == 0)
    {   
        int iup = (nprocess==1) ? sor->p_row + 1 : sor->p_row + 2;
        
        for(int i = 1; i < iup ; i++)
  	{
            for(int j = 0; j < sor->N; j++)
            {
                p_G[i][j] = G[i-1][j];
            }
        }
        for(int j = 0; j < sor->N; j++)
        {
            p_G[0][j] = 0.0; 
        }
        
        for(int k = 1; k < nprocess; k++)
        {
            int m_length = (k==nprocess-1) ? sor->rem_p_row + 1 : sor->p_row + 2;
            int begin    = (k * sor->p_row) - 1;
            
            MPI_Send(  G[begin], 
                        sor->N * m_length, 
                        MPI_DOUBLE, 
                        k, 
                        k, 
                        MPI_COMM_WORLD
                    );
        } 
    }
    else
    {
        int end = (rank==(nprocess-1)) ? sor->p_row+1 : sor->p_row+2;
        MPI_Recv(   p_G[0], 
                end * sor->N, 
                MPI_DOUBLE, 
                0, 
                rank, 
                MPI_COMM_WORLD, 
                MPI_STATUS_IGNORE
                );
    }
}

void sendSubMatricesFromSlavesToMaster( Sor* sor, const int rank, const int nprocess, 
                                        int p_G_lenght, int N, double p_G[p_G_lenght][N], 
                                        int M, double G[M][N]){

    if(rank == 0) 
    {
        int last = sor->p_row + 1; // same as p_G.length-1 -> (sor->p_row + 2) - 1 
        for(int i = 1; i < last; i++)
        {
            for(int j = 0; j < sor->N ; j++)
            {
                G[i-1][j] = p_G[i][j];
            }
        }
        int rm_length;
        for(int k = 1; k < nprocess; k++)
        {
            rm_length = (k == (nprocess-1)) ? sor->rem_p_row : sor->p_row;
            int first_line 	= k * sor->p_row;
            MPI_Recv( G[first_line], 
                        rm_length * sor->N, 
                        MPI_DOUBLE, 
                        k, 
                        k, 
                        MPI_COMM_WORLD, 
                        MPI_STATUS_IGNORE
                    );
        }
    } 
    else 
    {
        MPI_Ssend(  p_G[1], 
                    sor->N * sor->p_row, 
                    MPI_DOUBLE, 
                    0, 
                    rank, 
                    MPI_COMM_WORLD
                 );
    }
}

void swap_rows  (   int last_process, int master, int p_G_lenght, 
                    int N, double p_G[p_G_lenght][N], int rank)
{
           if(!last_process)
        {
            MPI_Sendrecv(   p_G[p_G_lenght-2], 
                            N,
                            MPI_DOUBLE, 
                            rank+1, 
                            1, 
                            p_G[p_G_lenght-1],  
                            N,
                            MPI_DOUBLE,
                            rank+1,
                            2, 
                            MPI_COMM_WORLD, 
                            MPI_STATUS_IGNORE
                        );
        }
        if(!master)
        {
            MPI_Sendrecv(   p_G[1], 
                            N,
                            MPI_DOUBLE, 
                            rank-1, 
                            2, 
                            p_G[0], 
                            N, 
                            MPI_DOUBLE,
                            rank-1,
                            1, 
                            MPI_COMM_WORLD, 
                            MPI_STATUS_IGNORE
                    );      
        } 
}

