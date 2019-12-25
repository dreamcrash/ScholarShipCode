/* 
 * File:   SeriesTest.h
 * Author: dreamcrash.
 *
 * Created on 2 de October de 2016, 16:46
 */

#ifndef SERIESTEST_H
#define SERIESTEST_H

void run                    (const int size, const int validation, int numThreads);
void JGFvalidate            (const int size, double TestArray[2][size]);
void Do                     (const int size, double TestArray [2][size]);
double TrapezoidIntegrate   (double x0, double x1, int nsteps, double omegan,  int select);
double thefunction          (double x, double omegan, int select);


#ifdef __cplusplus
extern "C" {
#endif




#ifdef __cplusplus
}
#endif

#endif /* SERIESTEST_H */

