/* 
 * File:   functions_mpi.h
 * Author: dreamcrash.
 *
 * Created on 13 de Outubro de 2016, 18:56
 */

#ifndef FUNCTIONS_MPI_H
#define FUNCTIONS_MPI_H

#include "Crypt.h"

int getProcessId        ();
int numberProcess       ();
void scattering         (Ideatest *ideatest, int rank, int numberProcesses);
void gathering          (Ideatest *ideatest, int rank, int numberProcesses);
#ifdef __cplusplus
extern "C" {
#endif




#ifdef __cplusplus
}
#endif

#endif /* FUNCTIONS_MPI_H */

