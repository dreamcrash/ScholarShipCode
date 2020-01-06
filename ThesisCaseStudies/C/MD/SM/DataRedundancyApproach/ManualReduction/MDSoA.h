/* 
 * File:   MDSoA.h
 * Author: Dreamcrash
 *
 */

#ifndef MDSOA_H
#define	MDSOA_H

#ifdef	__cplusplus
extern "C" {
#endif
#include "Structs.h"
#include "Reduction.h"

MD *newMD                       (int size,int *datasizes, int iterations, int numThreads);
void freeMD                     (MD *md);
void initialiseMD               (MD *md, int mm);
void printMDvariableControl     (MD *md);
void velocityXYZEscalar         (MD *md);
void cicleDoMove                (MD *md);
void cicleForces                (MD *md, Reduction **thrPrivateData);
void cicleMkekin                (MD *md);
void cicleVelavg                (MD *md);
void scaleTemperature           (MD *md,  int move);
void getFullPotentialEnergy     (MD *md,  int move);
void runMD                      (MD *md);

#ifdef	__cplusplus
}
#endif

#endif	/* MDSOA_H */

