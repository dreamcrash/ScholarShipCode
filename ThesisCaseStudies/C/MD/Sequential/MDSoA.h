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

MD *newMDSoA                       (int size,int *datasizes, int iterations);
void freeMDSoA                     (MD *md);
void initialiseMDSoA               (MD *md, int mm);
void printMDvariableControl        (MD *md);
void velocityXYZEscalarSoA         (MD *md);
void cicleDoMoveSoA                (MD *md);
void cicleForcesApproachSoAMPI     (MD *md);
void cicleMkekinSoA                (MD *md);
void cicleVelavgSoA                (MD *md);
void scaleTemperatureSoA           (MD *md,  int move);
void getFullPotentialEnergySoA     (MD *md,  int move);
void runMDSoA                      (MD *md);

#ifdef	__cplusplus
}
#endif

#endif	/* MDSOA_H */

