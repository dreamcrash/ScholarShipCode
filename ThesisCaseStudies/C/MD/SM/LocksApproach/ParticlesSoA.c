#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "Structs.h"
#include "ParticlesSoA.h"
#include "Random.h"

Particles *newParticles(int mdsize){
    
         size_t memoryNeeded = mdsize * sizeof(double);
    
         Particles *aux = malloc (sizeof(Particles));
             
         aux -> x  = malloc(memoryNeeded);
         aux -> y  = malloc(memoryNeeded);
         aux -> z  = malloc(memoryNeeded);
         aux -> vx = malloc(memoryNeeded);
         aux -> vy = malloc(memoryNeeded);
         aux -> vz = malloc(memoryNeeded);
         aux -> fx = malloc(memoryNeeded);
         aux -> fy = malloc(memoryNeeded);
         aux -> fz = malloc(memoryNeeded);
                        
         aux->numberParticles = mdsize;
         
         size_t memoryNeededLocking  = mdsize * sizeof(omp_lock_t);
         aux -> particlesLocks = malloc(memoryNeededLocking);   
         
         for(int i = 0; i < mdsize; i++)                       
         {
             omp_init_lock(&(aux -> particlesLocks[i]));
         }
         
         return aux;
}

void initParticles(Particles *p, int mdsize){
    
    for(int i = 0; i < mdsize; i++)
    {
        p->x[i]  = p->y[i]  = p->z[i]  = 0;
        p->vx[i] = p->vy[i] = p->vz[i] = 0;
        p->fx[i] = p->fy[i] = p->fz[i] = 0;                
    }
}

void freeParticles(Particles *p){
    
    free(p->fx);
    free(p->fy);
    free(p->fz);
    free(p->vx);
    free(p->vy);
    free(p->vz);
    free(p->x);
    free(p->y);
    free(p->z);
    
    for(int i = 0; i < p->numberParticles; i++)
    {
        omp_destroy_lock(&(p -> particlesLocks[i]));
    }
    free(p->particlesLocks);
}

/** Scales the velocity of all particles **/
void scaleVelocity(Particles *particles, int numParticles, double sc){
    
    for(int i = 0; i < numParticles; ++i)
    {
        particles->vx[i] *= sc;
        particles->vy[i] *= sc;
        particles->vz[i] *= sc;
    }
}

double getAverageVelocity(Particles *particle){
    
    int numberParticles = particle->numberParticles;
    double vel = 0.0;
    
    for (int pos = 0; pos < numberParticles; pos++)
    {
        vel += sqrt(
                        particle->vx[pos] * particle->vx[pos] +
                        particle->vy[pos] * particle->vy[pos] +
                        particle->vz[pos] * particle->vz[pos] 
                    );
    }
    return vel;
}
    
void printParticles(Particles *p, int numParticles){
    
    for (int i = 0; i < numParticles; i++)
    {
        printf("------ Particle %d ------\n",i);
        printf("P(x,y,z) = (%f,%f,%f)\n",p->x[i] ,p->y[i] ,p->z[i]);
	printf("V(x,y,z) = (%f,%f,%f)\n",p->vx[i],p->vy[i],p->vz[i]);
	printf("F(x,y,z) = (%f,%f,%f)\n",p->fx[i],p->fy[i],p->fz[i]);
    }
    
}

/** Does the movement for all particles in the system. */
void calculate_move(Particles *p, double side){

    int numberParticles = p->numberParticles;
    for (int i = 0; i < numberParticles; i++)
    {
	domove(p, side, i);
    }
}

/** Moving one particle */
void domove(Particles *p, double side, int pos){
    
     p->x[pos] += p->vx[pos] + p->fx[pos];
     p->y[pos] += p->vy[pos] + p->fy[pos];
     p->z[pos] += p->vz[pos] + p->fz[pos];

     if 	(p->x[pos] < 0)    p->x[pos] += side;
     else if    (p->x[pos] > side) p->x[pos] -= side;
     if 	(p->y[pos] < 0)    p->y[pos] += side;
     else if    (p->y[pos] > side) p->y[pos] -= side;
     if 	(p->z[pos] < 0)    p->z[pos] += side;
     else if    (p->z[pos] > side) p->z[pos] -= side;

     p->vx[pos] += p->fx[pos];
     p->vy[pos] += p->fy[pos];
     p->vz[pos] += p->fz[pos];
		
     p->fx[pos] = p->fy[pos] = p->fz[pos] = 0.0;
}

void calculate_force(MD *md){
    
    int numberParticles = md->mdsize;
    #pragma omp parallel for schedule(static, 1)
    for(int i = 0; i < numberParticles; i ++)
    {
        force3Law(md, md->particlesSOA, i);
    }
}

double mkekin(Particles *p, double hsq2){
    
    int numberParticles = p->numberParticles;
    double sum = 0.0;
    
    for (int pos = 0; pos < numberParticles; pos++)
    {
        p->fx[pos] *= hsq2;
        p->fy[pos] *= hsq2;
        p->fz[pos] *= hsq2;
        
        p->vx[pos] += p->fx[pos];
        p->vy[pos] += p->fy[pos];
        p->vz[pos] += p->fz[pos];
        
        sum +=  p->vx[pos] * p->vx[pos] + 
                p->vy[pos] * p->vy[pos] + 
                p->vz[pos] * p->vz[pos];
    }
    
    return sum;
}

/** Generating particles position **/
void particleGenerate(Particles *p, int mm, double a){
    
        int lg,i,k,j,particle = 0;
        for (lg = 0; lg <= 1; lg++) 
         for (i = 0; i < mm; i++) 
          for (j = 0; j < mm; j++) 
           for (k = 0; k < mm; k++)
           {
                p->x[particle] = i * a + lg * a * 0.5;
	        p->y[particle] = j * a + lg * a * 0.5; 
	        p->z[particle] = k * a;				 
	        particle++;
           }
        for (lg = 1; lg <= 2; lg++) 
	 for (i = 0; i < mm; i++) 
          for (j = 0; j < mm; j++) 
	   for (k = 0; k < mm; k++)
           {
                p->x[particle] = i * a + (2 - lg) * a * 0.5;
		p->y[particle] = j * a + (lg - 1) * a * 0.5;
                p->z[particle] = k * a + a * 0.5;
	        particle++;
	   }
}	

/** Generating Particles Velocity */
void velocityGenerate(Particles *p, int numOfParticles){
    
      double r = 0.0;
      
      Random *randnum = newRandom(0, 0.0,0.0);
		
      for (int i = 0; i < numOfParticles; i += 2) 
      {
	  r = seed(randnum);
	  p->vx[i]      = r * randnum->v1;
	  p->vx[i + 1]  = r * randnum->v2;
      }

      for (int i = 0; i < numOfParticles; i += 2)
      {
	  r = seed(randnum);
	  p->vy[i]      = r * randnum->v1;
	  p->vy[i + 1]  = r * randnum->v2;
      }

      for (int i = 0; i < numOfParticles; i += 2)
      {
	  r = seed(randnum);
	  p->vz[i]      = r * randnum->v1;
	  p->vz[i + 1]  = r * randnum->v2;
      }
      
      free(randnum);
  }

/** Scaling Velocity component */
double scalingVelocity(Particles *p, int numOfParticles){
    
     double sp = 0.0, ekin = 0.0;
     
     /** Scaling vx */
     for(int i = 0; i < numOfParticles; i++)	
     {
	sp += p->vx[i];
     }
	    
     sp /= numOfParticles;
	    
     for (int i = 0; i < numOfParticles; i++) 
     {
        p->vx[i] -= sp;
	ekin += p->vx[i] * p->vx[i];
     }
     
     sp = 0.0;
     
     /** Scaling vy */
     for(int i = 0; i < numOfParticles; i++)	
     {
	sp += p->vy[i];
     }
	    
     sp /= numOfParticles;
	    
     for (int i = 0; i < numOfParticles; i++) 
     {
        p->vy[i] -= sp;
	ekin += p->vy [i]* p->vy[i];
     }
     
     sp = 0.0;
     
     /** Scaling vz*/
     for(int i = 0; i < numOfParticles; i++)	
     {
	sp += p->vz[i];
     }
	    
     sp /= numOfParticles;
	    
     for (int i = 0; i < numOfParticles; i++) 
     {
        p->vz[i] -= sp;
	ekin += p->vz[i] * p->vz[i];
     }
     
     return ekin;
}
	
/** Calculating the force of one particle with the remaining */
void force3Law(MD *md, Particles *p, int particleA){
    
        const double radius = md->rcoffs;
        const double side = md->side;
        const double sideh = 0.5  * md->side;
        const int numParticles = md->mdsize;
        
        /** Get particle A position **/
	const double xi = p->x[particleA];		
	const double yi = p->y[particleA];
	const double zi = p->z[particleA];
        
        /** Accumulate Particle A force **/
	double fxAcc = 0, fyAcc = 0, fzAcc = 0;
	
	/** Save MD temporary variables values **/
	double epot = 0.0, vir = 0.0;
	int interactions = 0;
        
	/** Calculate the force interaction between pairs of particles */
        for (int particleB = particleA + 1; particleB < numParticles; particleB++)
	{
            /** Calculate the distance between particle A and particle particleB **/
            double xx = xi - p->x[particleB];
            double yy = yi - p->y[particleB];
            double zz = zi - p->z[particleB];
            
            /** Check if the distance is inside the world **/
            if      (xx < (-sideh)) 	{ xx += side; }
	    else if (xx > (sideh))  	{ xx -= side; }
	    if      (yy < (-sideh)) 	{ yy += side; }
	    else if (yy > (sideh))  	{ yy -= side; }
	    if      (zz < (-sideh)) 	{ zz += side; }
	    else if (zz > (sideh))  	{ zz -= side; }

            const double rd = xx*xx + yy*yy + zz*zz;
	           
	    if(rd <= radius)
	    {
                const double rrd = 1.0/rd;
	        const double rrd2 = rrd  * rrd;
	        const double rrd3 = rrd2 * rrd;
	        const double rrd4 = rrd2 * rrd2;
	        const double rrd6 = rrd2 * rrd4;
	        const double rrd7 = rrd6 * rrd;
	        const double r148 = rrd7 - 0.5 * rrd4;
	         
                const double tmpFx = xx * r148;       
                const double tmpFy = yy * r148;
                const double tmpFz = zz * r148;
                 
                /** Force Accumulation fx,fy,fz*/
                fxAcc += tmpFx;
                fyAcc += tmpFy;
                fzAcc += tmpFz;
                 
                /** Calculating thirds Newton's law */
                omp_set_lock(&(p->particlesLocks[particleB]));     
                p->fx[particleB] -= tmpFx;
                p->fy[particleB] -= tmpFy;
                p->fz[particleB] -= tmpFz;  
                omp_unset_lock(&(p->particlesLocks[particleB]));
                
    		epot += rrd6 - rrd3;
		vir  -= rd * r148;
		interactions++;             
	    }
        }
        /** Force Calculation **/
        omp_set_lock(&(p->particlesLocks[particleA]));
        p->fx[particleA] += fxAcc;
        p->fy[particleA] += fyAcc;
        p->fz[particleA] += fzAcc;
        omp_unset_lock(&(p->particlesLocks[particleA]));
        
        /** MD control variables actualization */
        #pragma omp critical(criticalMDvariables)
        {
            md->epot         += epot;
            md->vir          += vir;
            md->interactions += interactions;
        }
}


