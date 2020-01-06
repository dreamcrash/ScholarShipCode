/* 
 * File:   Structs.h
 * Author: Dreamcrash
 *
 */

#ifndef STRUCTS_H
#define	STRUCTS_H

#ifdef	__cplusplus
extern "C" {
#endif

// Using a SoA Layout   
typedef struct PARTICLES{

    double *x,  *y,  *z;  // Positions
    double *vx, *vy, *vz; // Velocities
    double *fx, *fy, *fz; // Forces
    int numberParticles;    
    
} Particles;

typedef struct MD{
    
    Particles *particlesSOA;   // Use for the layouts SoA

    int irep;
    int istop;    
    int iprint;
    int mdsize;
    int size;
    int movemx;
    int interactions;
    double den;
    double tref;
    double h;
    double epot, vir;
    double rcoffs, side, hsq, hsq2, vel;
    double sum, tscale, sc, ekin, ek;
    
    //DM CCCs
    int processID;
    int totalProcess;
    
    // SM CCCs
    int numThreads;
}MD;

#ifdef	__cplusplus
}
#endif

#endif	/* STRUCTS_H */

