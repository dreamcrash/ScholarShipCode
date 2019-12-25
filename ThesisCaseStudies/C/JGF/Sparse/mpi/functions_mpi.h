/* 
 * File:   functions_mpi.h
 * Author: Dreamcrash
 *
 * Created on 13 de October de 2016, 18:56
 */

#ifndef FUNCTIONS_MPI_H
#define FUNCTIONS_MPI_H

int getProcessId            ();
int numberProcess           ();
int master                  ();
void send_data_to_slaves    (int *row, int *col, double *val, 
                                int dim, int slaveID);
void recv_data_from_master  (int *row, int *col, double *val, 
                                int dim, int slaveID);


#ifdef __cplusplus
extern "C" {
#endif




#ifdef __cplusplus
}
#endif

#endif /* FUNCTIONS_MPI_H */

