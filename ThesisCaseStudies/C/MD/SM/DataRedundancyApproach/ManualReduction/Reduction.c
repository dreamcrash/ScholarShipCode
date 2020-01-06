#include <stdio.h>
#include <stdlib.h>
#include "Structs.h"
#include "Reduction.h"

/** Create structure to reduce variables **/
Reduction **createReduction(Particles *p, int numOfParticles, int numThreads){
   
    size_t memory = numOfParticles * sizeof(double);
    Reduction **r = malloc(numThreads * sizeof(Reduction*));
    
    /** Master **/
    r[0] = malloc(sizeof(Reduction));
    r[0]->epot = 0.0;
    r[0]->vir  = 0.0;
    r[0]->interactions = 0;
    r[0]->fx = p->fx; // using the original data
    r[0]->fy = p->fy; // using the original data
    r[0]->fz = p->fz; // using the original data  
        
    /** Slaves */
    for(int i = 1; i < numThreads; i++)
    {
        r[i] = malloc(sizeof(Reduction));
        r[i]->epot = 0.0;
        r[i]->vir  = 0.0;
        r[i]->interactions = 0;
        r[i]->fx = malloc(memory);
        r[i]->fy = malloc(memory);
        r[i]->fz = malloc(memory);
    }
    return r;
}

/** Free Reduction Variables Memory **/
void freeReductionVariables(Reduction **vars, int numThreads){
    
    free(vars[0]);  // Free Master
    
    /** Free the slaves data */
    for(int i = 1; i < numThreads; i++)
    {
        free(vars[i]->fx);
        free(vars[i]->fy);
        free(vars[i]->fz);
        free(vars[i]);
    }
    free(vars);
}

/** Initialization of Variable Reductions */
void initForces(Reduction *vars, int numOfParticles){
    
    for(int i = 0; i < numOfParticles; i++)
    {
        vars->fx[i] = vars->fy[i] = vars->fz[i] = 0;
    }
}

/** Initialization of Variable Reductions */
void initMDVariables(Reduction *vars){
    
    vars->epot = vars->vir = 0.0;
    vars->interactions = 0;
}

void initThreadRelatedVariables(Reduction **vars, int threadID, int numberParticles){
    
    initMDVariables(vars[threadID]);
    if(threadID != 0)
    {
        initForces(vars[threadID], numberParticles);
    }
}

void reduceThreadRelatedVariables(MD *md, Reduction **vars) {

    md->epot += vars[0]->epot;
    md->vir += vars[0]->vir;
    md->interactions += vars[0]->interactions;

    for (int t = 1; t < md->numThreads; t++) 
    {
        for (int p = 0; p < md->particlesSOA->numberParticles; p++) 
        {
            md->particlesSOA->fx[p] += vars[t]->fx[p];
            md->particlesSOA->fy[p] += vars[t]->fy[p];
            md->particlesSOA->fz[p] += vars[t]->fz[p];
        }
        md->epot += vars[t]->epot;
        md->vir += vars[t]->vir;
        md->interactions += vars[t]->interactions;
    }  
}