/* 
 * File:   RandomMatrix.c
 * Author: Dreamcrash
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "RandomMatrix.h"


/**
 * Fill up matrix with random values
 * @param min   : The lowest value that can be generate
 * @param max   : The highest value that can be generate
 * @param seed  : The seed of the generator
 * @param rows   : The total of rows
 * @param cols   : The total of columns
 * @param m     : The matrix to be fill
 */
void fillupRandomly (const int min, const int max, const int seed, 
                         const int rows, const int cols, double m[rows][cols]){
     srand(seed);

    for (int i = 0; i < rows; i++)
        for (int j = 0; j < cols; j++)
            m[i][j] = (((double) rand() / (double) RAND_MAX)) * (max - min) + min;
} 
