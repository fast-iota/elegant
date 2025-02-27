/*************************************************************************\
* Copyright (c) 2017 The University of Chicago, as Operator of Argonne
* National Laboratory.
* This file is distributed subject to a Software License Agreement found
* in the file LICENSE that is included with this distribution. 
\*************************************************************************/

/*
 * Joe Calvey, Michael Borland 2017
 */

typedef struct {
  long nGasses;
  char **gasName;
  long nLocations;
  double *s;         /* s[j] is the location of the jth set of pressure samples */
  double **pressure; /* pressure[i][j] is the pressure of the ith species at the jth location */
  double temperature; /* in degrees K */
} PRESSURE_DATA;

#ifdef __cplusplus
extern "C" {
#endif

void readGasPressureData(char *filename, PRESSURE_DATA *pressureData, double factor);
void computeAverageGasPressures(double sStart, double sEnd, double *pressure, PRESSURE_DATA *pressureData);

#ifdef __cplusplus
}
#endif
