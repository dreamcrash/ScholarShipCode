/* 
 * File:   functions_mpi.h
 * Author: Dreamcrash
 *
 */

#ifndef FUNCTIONS_MPI_H
#define	FUNCTIONS_MPI_H

#ifdef	__cplusplus
extern "C" {
#endif
#include "Structs.h"
    
int getProcessId                        ();
int numberProcess                       ();
void reduceForces                       (MD *md);
void reduceStaticVariable               (MD *md);



#ifdef	__cplusplus
}
#endif

#endif	/* FUNCTIONS_MPI_H */

