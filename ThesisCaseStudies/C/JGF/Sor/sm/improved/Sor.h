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

#define CACHELINE 128

void run            (int size, int validation, int nthreads);
void JGFinitialise  (Sor *sor);
void JGFKernel      (Sor *sor, int total_threads);
void RandomMatrix   (int M, int N, double G [M][N]);
void sor_simulation (double omega, int M, int N, double G[M][N], int iterations, 
                        int total_threads, volatile long sync[total_threads][CACHELINE]);
void row_iterations (int begin, int end, double omega_over_four, double one_minus_omega, 
                        int Mm1, int Nm1, double G[Mm1+1][Nm1+1]);
void call_barrier(int id_thread, int total_threads, volatile long sync[total_threads][CACHELINE]);

void middle_rows    (const double omega_over_four, const double one_minus_omega, 
                        const int Nm1, double *Gi, double *Gim1, 
                        double *G1, double *G2); 
void last_row       (const double omega_over_four, const double one_minus_omega, 
                     const int Nm1, double *Gi, double *Gim1, double *Gim2_);

void first_row      ( const double omega_over_four, const double one_minus_omega, 
                        const int Nm1, double *Gi, double *Gim1, double *Gip1_);


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

