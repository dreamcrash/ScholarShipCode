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

Particles *newParticles          (int mdsize);
void initParticles               (Particles *p, int mdsize);
void freeParticles               (Particles *p);
void scaleVelocity               (Particles *p, int numParticles, double sc);
double getAverageVelocity        (Particles *particle);
void printParticles              (Particles *p, int numParticles);
void calculate_move              (Particles *p, double side);
void domove                      (Particles *p , double side, int pos);
void calculate_force             (MD *md);
double mkekin                    (Particles *p, double hsq2);
void particleGenerate            (Particles *p, int mm, double a);
void velocityGenerate            (Particles *p, int numOfParticles);
double scalingVelocity           (Particles *p, int numOfParticles);
void force3Law                   (MD *md, Particles *p, int particleA);


#ifdef	__cplusplus
}
#endif

#endif	/* PARTICLESSOA_H */

