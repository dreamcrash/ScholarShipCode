#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "Structs.h"
#include "MDSoA.h"
#include "ParticlesSoA.h"

MD *newMD(int size,int *datasizes, int iterations){
    
    int mm;
    
    MD *md = malloc(sizeof (MD));
    md->size = size;
    md->movemx = iterations;
    mm = datasizes[md->size];
    md->mdsize = mm * mm * mm * 4;
    
    /** Creating the Particles structure */
    md->particlesSOA = newParticles(md->mdsize);       
    initialiseMD(md, mm);

    return md;
}

void freeMD(MD *md){
    
    freeParticles(md->particlesSOA);
    free(md->particlesSOA);
    free(md);
}

void printMDvariableControl(MD *md){
    
    printf("--------- MD ---------\n");
    printf("mdsize = %d\n",md->mdsize);
    printf("ek = %f\n",md->ek);
    printf("ekin = %f\n",md->ekin);
    printf("sc = %f\n",md->sc);
    printf("sum = %f\n",md->sum);
    printf("vel = %f\n",md->vel);
    printf("epot = %f\n",md->epot);
    printf("vir = %f\n",md->vir);
    printf("interactions = %d\n",md->interactions);
}

void initialiseMD(MD *md, int mm){
    
    double a, rcoff;
   
    md->irep    = 10;
    md->istop   = 19;
    md->iprint  = 10;
    md->den     = 0.83134;
    md->tref    = 0.722;
    md->h       = 0.064;
    md->epot    = 0.0;
    md->vir     = 0.0;
    md->vel     = 0.0;
    md->interactions = 0;
    
    md->side = pow((md->mdsize / md->den), 0.3333333);
    rcoff = mm / 4.0;
    a = md->side / mm;
    md->hsq = md->h * md->h;
    md->hsq2 = md->hsq * 0.5;
    md->rcoffs = rcoff * rcoff;
    md->tscale = 16.0 / (1.0 * md->mdsize - 1.0);
    
    initParticles       (md->particlesSOA, md->mdsize);
    particleGenerate    (md->particlesSOA, mm, a);
    velocityGenerate    (md->particlesSOA, md->mdsize);
    velocityXYZEscalar  (md);
}

/** Scaling Particles Velocities*/
void velocityXYZEscalar(MD *md){
    
    double ts;
    md->ekin = 0.0;
    
    md->ekin += scalingVelocity (md->particlesSOA,md->mdsize);
    
    ts = md->tscale * md->ekin;
    md->sc = md->h * sqrt(md->tref / ts);
    
    scaleVelocity(md->particlesSOA,md->mdsize,md->sc);
}

/** Move the particles*/
void cicleDoMove(MD *md){
    calculate_move(md->particlesSOA, md->side);
}

/** Compute forces with 3 Newton's Law*/ 
void cicleForcesApproach(MD *md){

    md->epot = md->vir = 0.0;
    
    calculate_force (md);
}

/** scale forces, update velocities */
void cicleMkekin(MD *md){
    md->sum = mkekin(md->particlesSOA, md->hsq2);
}

/** Average velocity */
void cicleVelavg(MD *md){
    
    md->ekin = md->sum / md->hsq;
    md->vel = 0.0;
	     
    md->vel = getAverageVelocity(md->particlesSOA);
    
    md->vel /= md->h;
}

/** Scaling temperature */	
void scaleTemperature(MD *md,  int move){
    
    if ((move < md->istop) && (((move + 1) % md->irep) == 0))
    {
       md->sc = sqrt(md->tref / (md->tscale * md->ekin));
       scaleVelocity(md->particlesSOA, md->mdsize, md->sc);  // Scaling velocities
       md->ekin = md->tref / md->tscale;
    }
}	

/** Getting full potential energy */
void getFullPotentialEnergy(MD *md, int move){
    
    if (((move + 1) % md->iprint) == 0) 
    {
        md->ek   = 24.0 * md->ekin;
	md->epot = 4.0  * md->epot;
	md->vel /= md->mdsize;
    }
}

/** MD Simulation */
void runMD(MD *md){
     
    for(int move = 0; move < md->movemx; ++move)
    {
        cicleDoMove                  (md);
        cicleForcesApproach          (md);
        cicleMkekin                  (md);
        cicleVelavg                  (md);
        scaleTemperature             (md,move);
        getFullPotentialEnergy       (md,move);
    }   
}
