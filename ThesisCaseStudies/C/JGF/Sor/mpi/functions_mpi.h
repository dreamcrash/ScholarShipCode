/* 
 * File:   functions_mpi.h
 * Author: Dreamcrash
 *
 * Created on 13 of October of 2016, 18:56
 */

#ifndef FUNCTIONS_MPI_H
#define FUNCTIONS_MPI_H
#include "Sor.h"

int getProcessId                        ();
int numberProcess                       ();
void sendMatrixFromMasterToSlaves       (  Sor* sor, const int rank, const int nprocess, 
                                            int p_G_lenght, int N, double p_G[p_G_lenght][N], 
                                            int M, double G[M][N]);
void sendSubMatricesFromSlavesToMaster  ( Sor* sor, const int rank, const int nprocess, 
                                            int p_G_lenght, int N, double p_G[p_G_lenght][N], 
                                            int M, double G[M][N]);
void swap_rows  (   int last_process, int master, int p_G_lenght, 
                    int N, double p_G[p_G_lenght][N], int rank);

#ifdef __cplusplus
extern "C" {
#endif




#ifdef __cplusplus
}
#endif

#endif /* FUNCTIONS_MPI_H */

