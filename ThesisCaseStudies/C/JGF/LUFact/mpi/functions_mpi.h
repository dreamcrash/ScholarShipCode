/* 
 * File:   functions_mpi.h
 * Author: Dreamcrash.
 *
 * Created on 13 of October of 2016, 18:56
 */

#ifndef FUNCTIONS_MPI_H
#define FUNCTIONS_MPI_H

int getProcessId                        ();
int numberProcess                       ();
void sendAndRecvSubMatrixCycleManner    (int max_row, int max_col, double **a);
void send_buffer_to_worker              (double *buf_col_k, int n, int worker);
void send_l_to_worker                   (int *tmp_l, int worker);
void master_collecting_results          (int totalProcess, int max_row, int max_col, double **a);
void slaves_send_results                (int slaveID, int totalProcess, int max_row, int max_col, double **a);


#ifdef __cplusplus
extern "C" {
#endif




#ifdef __cplusplus
}
#endif

#endif /* FUNCTIONS_MPI_H */

