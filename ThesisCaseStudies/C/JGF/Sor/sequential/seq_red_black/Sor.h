/* 
 * File:   Sor.h
 * Author: Dreamcrash
 *
 * Created on 20 of October of 2016, 15:58
 */

#ifndef SOR_H
#define SOR_H


typedef struct SOR{
    
    int     JACOBI_NUM_ITER;
    long    RANDOM_SEED;
    
    int     size;
    
    int M, N;
    
    double Gtotal;
    double ops;
    
}Sor;

void run            (int size, int validation);
void JGFinitialise  (Sor *sor);
void JGFKernel      (Sor *sor);
void RandomMatrix   (int M, int N, double G [M][N]);
void sor_simulation (double omega, int M, int N, double G[M][N], int iterations);

int abs_            (int value);
double abs_double   (double value);
void JGFvalidate    (double Gtotal, int size);


#ifdef __cplusplus
extern "C" {
#endif


#ifdef __cplusplus
}
#endif

#endif /* SOR_H */

