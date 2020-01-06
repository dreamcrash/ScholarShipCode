#include "mpi.h"
#include "Structs.h"
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

void reduceForces(MD *md){
    
    int size = md->mdsize;
    MPI_Allreduce(MPI_IN_PLACE, md->particlesSOA->fx, size,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, md->particlesSOA->fy, size,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, md->particlesSOA->fz, size,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
}

void reduceStaticVariable(MD *md){
    
    MPI_Allreduce(MPI_IN_PLACE	,&(md->epot)           ,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE	,&(md->vir)            ,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE	,&(md->interactions)   ,1,MPI_INT   ,MPI_SUM,MPI_COMM_WORLD);
}



