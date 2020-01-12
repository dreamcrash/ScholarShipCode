/* 
 * File:   main.h
 * Author: Dreamcrash
 *
 */

#ifndef MAIN_H
#define MAIN_H

#ifdef __cplusplus
extern "C" {
#endif

void validateMM(    const int maxRowA, const int maxColB, const int maxColA,
                    double A[maxRowA][maxColA], 
                    double B[maxColA][maxColB], 
                    double C[maxRowA][maxColB]);


#ifdef __cplusplus
}
#endif

#endif /* MAIN_H */

