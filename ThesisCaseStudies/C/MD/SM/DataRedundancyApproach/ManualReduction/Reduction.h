/* 
 * File:   Reduction.h
 * Author: Dreamcrash
 *
 */

#ifndef REDUCTION_H
#define	REDUCTION_H

#ifdef	__cplusplus
extern "C" {
#endif

typedef struct REDUCTION{
    
    double epot;
    double vir;
    int interactions;
    double *fx, *fy, *fz;
    
}Reduction;    

typedef struct MD_R{
    
    double epot;
    double vir;
    int interactions;
    
}Md_r;    
    
Reduction **createReduction         (Particles *p, int numOfParticles, int numThreads);
void freeReductionVariables         (Reduction **thrPrivateData, int numThreads);
void reductionSoA                   (Reduction **thrPrivateData, MD *md);
void initForces                     (Reduction *thrPrivateData, int numOfParticles);
void initMDVariables                (Reduction *thrPrivateData);
void initThreadRelatedVariables     (Reduction **thrPrivateData, int threadID, int numOfParticles);
void reduceThreadRelatedVariables   (MD *md, Reduction **vars);


#ifdef	__cplusplus
}
#endif

#endif	/* REDUCTION_H */