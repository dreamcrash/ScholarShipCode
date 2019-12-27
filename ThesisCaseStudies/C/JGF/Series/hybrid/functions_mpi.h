/* 
 * File:   functions_mpi.h
 * Author: Dreamcrash
 *
 * Created on 2 de October de 2016, 17:08
 */

#ifndef FUNCTIONS_MPI_H
#define FUNCTIONS_MPI_H

int getProcessId        ();
int numberProcess       ();
int master              ();
int slaves              ();
void reduceTestArray    (   const int process, const int size, 
                            double TestArray [2][size]);

#ifdef __cplusplus
extern "C" {
#endif




#ifdef __cplusplus
}
#endif

#endif /* FUNCTIONS_MPI_H */

