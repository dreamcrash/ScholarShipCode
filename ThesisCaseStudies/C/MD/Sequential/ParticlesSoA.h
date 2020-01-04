/* 
 * File:   ParticlesSoA.h
 * Author: Dreamcrash
 *
 */

#ifndef PARTICLESSOA_H
#define	PARTICLESSOA_H

#ifdef	__cplusplus
extern "C" {
#endif
#include "Structs.h"

Particles *newParticlesSoA          (int mdsize);
void initParticlesSoA               (Particles *p, int mdsize);
void freeParticlesArraysSoA         (Particles *p);
void scaleVelocitySoA               (Particles *p, int numParticles, double sc);
double getAverageVelocity           (Particles *particle);
void printParticlesSoA              (Particles *p, int numParticles);
void calculate_move                 (Particles *p, double side);
void domoveSoA                      (Particles *p , double side, int pos);
void calculate_force                (MD *md);
double mkekin                       (Particles *p, double hsq2);
void particleGenerateSoA            (Particles *p, int mm, double a);
void velocityGenerateSoA            (Particles *p, int numOfParticles);
double scalingVelocitySoA           (Particles *p, int numOfParticles);
void force3LawSoAReduction          (MD *md, int particleA, Particles *p);


#ifdef	__cplusplus
}
#endif

#endif	/* PARTICLESSOA_H */

