/* 
 * File:   processMPI.h
 * Author: Dreamcrash
 *
 */

#ifndef PROCESSMPI_H
#define PROCESSMPI_H

#ifdef __cplusplus
extern "C" {
#endif

int  getProcessId   ();
int  numberProcess  ();
int  myMatrixSize (const int totalSize, const int processID, const int numProcesses);
void masterSendSubMatricesToSlaves (  const int maxRow, const int maxCol, double m[][maxCol], 
                                        const int tile, const int idProcess, const int numberProcesses);
void masterCollectsSubMatrix( const int maxRow, const int maxCol, double m[][maxCol], 
                                        const int tile, const int idProcess, const int numberProcesses);


#ifdef __cplusplus
}
#endif

#endif /* PROCESSMPI_H */

