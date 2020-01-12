/* 
 * File:   RandomMatrix.h
 * Author: Dreamcrash
 *
 */

#ifndef RANDOMMATRIX_H
#define RANDOMMATRIX_H

#ifdef __cplusplus
extern "C" {
#endif
 
#define MAX 5
#define MIN -5
#define RANDOM_SEED 123456

void fillupRandomly (   const int min, const int max, const int seed, 
                        const int rows, const int cols, double m[rows][cols]);

#ifdef __cplusplus
}
#endif

#endif /* RANDOMMATRIX_H */

