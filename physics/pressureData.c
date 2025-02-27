/*************************************************************************\
* Copyright (c) 2017 The University of Chicago, as Operator of Argonne
* National Laboratory.
* This file is distributed subject to a Software License Agreement found
* in the file LICENSE that is included with this distribution. 
\*************************************************************************/

/* file: ionEffects.c
 * purpose: simulation of ion interaction with the beam
 *
 * Joe Calvey, Michael Borland 2017
 */
#if defined(SOLARIS) && !defined(__GNUC__)
#include <sunmath.h>
#endif

#include "mdb.h"
#include "SDDS.h"
#include "constants.h"
#include "pressureData.h"

long checkSddsColumn(SDDS_TABLE *SDDS_table, char *name, char *units)
{
  char *units1;
  if (SDDS_GetColumnIndex(SDDS_table, name)<0)
    return(0);
  if (SDDS_GetColumnInformation(SDDS_table, "units", &units1, SDDS_GET_BY_NAME, name)!=SDDS_STRING) {
    SDDS_SetError("units field of column has wrong data type!");
    SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors|SDDS_VERBOSE_PrintErrors);
  }
  if (!units || SDDS_StringIsBlank(units)) {
    if (!units1)
      return(1);
    if (SDDS_StringIsBlank(units1)) {
      free(units1);
      return(1);
    }
    return(0);
  }
  if (!units1)
    return(0);
  if (strcmp(units, units1)==0) {
    free(units1);
    return(1);
  }
  free(units1);
  return(0);
}

void readGasPressureData(char *filename, PRESSURE_DATA *pressureData)
{
  /* Assumed file structure:
   * Parameters: 
   * Gasses --- SDDS_STRING giving comma- or space-separated list of gas species, e.g., "H2O H2 N2 O2 CO2 CO CH4"
   * Temperature --- SDDS_FLOAT or SDDS_DOUBLE giving temperature in degrees K. Defaults to 293.
   * Columns:
   * s         --- SDDS_FLOAT or SDDS_DOUBLE giving location in the lattice
   * <gasName> --- SDDS_FLOAT or SDDS_DOUBLE giving pressure of <gasName> in Torr or nT
   */

  SDDS_DATASET SDDSin;
  char *gasColumnList, *ptr;
  long i;
  double dsMin, dsMax, ds, pressureMultiplier;

  if (!SDDS_InitializeInput(&SDDSin, filename))
    bombVA("Failed to initialize input from %s", filename);

  if (!checkSddsColumn(&SDDSin, "s", "m"))
    bombVA("Column 's' is missing or does not have units of 'm' in %s", filename);
  if (SDDS_CheckParameter(&SDDSin, "Gasses", NULL, SDDS_STRING, stdout)!=SDDS_CHECK_OK)
    bombVA("Parameters \"Gasses\" is missing or not string type in %s", filename);

  if (SDDS_ReadPage(&SDDSin)<=0) 
    SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors|SDDS_VERBOSE_PrintErrors);

  if (!SDDS_GetParameter(&SDDSin, "Gasses", &gasColumnList))
    SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors|SDDS_VERBOSE_PrintErrors);
    
  switch (SDDS_CheckParameter(&SDDSin, "Temperature", "K", SDDS_ANY_FLOATING_TYPE, NULL)) {
  case SDDS_CHECK_OK:
    if (!SDDS_GetParameterAsDouble(&SDDSin, "Temperature", &pressureData->temperature) ||
        pressureData->temperature<=0) 
      bombVA("Problem reading 'Temperature' from %s. Check for valid value (got %le).\n", filename,
                    pressureData->temperature);
    break;
  case SDDS_CHECK_NONEXISTENT:
    pressureData->temperature = 273+20;
    printf("Parameter 'Temperature' missing from %s, assuming %le K\n", filename, pressureData->temperature);
    break;
  case SDDS_CHECK_WRONGTYPE:
    bombVA("Parameter 'Temperature' in %s has wrong type. Expect SDDS_DOUBLE or SDDS_FLOAT.\n", filename);
    break;
  case SDDS_CHECK_WRONGUNITS:
    bombVA("Parameter 'Temperature' in %s has wrong units. Expect 'K'.\n", filename);
    break;
  default:
    bombVA("Unexpected value checking existence, units, and type for 'Temperature' in %s\n", filename);
    break;
  }
  
  pressureData->nGasses = 0;
  pressureData->gasName = NULL;
  while ((ptr=get_token(gasColumnList))!=NULL) {
    pressureData->gasName = (char**)SDDS_Realloc(pressureData->gasName, sizeof(*(pressureData->gasName))*(pressureData->nGasses+1));
    cp_str(&pressureData->gasName[pressureData->nGasses], ptr);
    pressureData->nGasses += 1;
  }
  free(gasColumnList);
  pressureData->gasData = tmalloc(sizeof(*(pressureData->gasData))*pressureData->nGasses);
  for (i=0; i<pressureData->nGasses; i++) {
    if (!identifyGas(pressureData->gasData+i, pressureData->gasName[i])) {
      fprintf(stderr, "unknown gas: \"%s\"\n", pressureData->gasName[i]);
    }
  }

  pressureData->nLocations = SDDS_RowCount(&SDDSin);
  /* printf("Gas data provided at %ld s locations\n", pressureData->nLocations); */
  if (!(pressureData->s=SDDS_GetColumnInDoubles(&SDDSin, "s")))
    SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors|SDDS_VERBOSE_PrintErrors);

  pressureData->pressure = (double**)czarray_2d(sizeof(double), pressureData->nLocations, pressureData->nGasses);
  for (i=0; i<pressureData->nGasses; i++) {
    pressureMultiplier = 1;
    if (!checkSddsColumn(&SDDSin, pressureData->gasName[i], "Torr") && !checkSddsColumn(&SDDSin, pressureData->gasName[i], "T")) {
      pressureMultiplier = 1e-9;
      if (!checkSddsColumn(&SDDSin, pressureData->gasName[i], "nT") && !checkSddsColumn(&SDDSin, pressureData->gasName[i], "nTorr"))
        bombVA("Column \"%s\" is missing, not floating-point type, or does not have units of \"Torr\" or \"nT\" in %s", 
                      pressureData->gasName[i], filename);
    }
    if (!(pressureData->pressure[i] = SDDS_GetColumnInDoubles(&SDDSin, pressureData->gasName[i]))) {
      SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors|SDDS_VERBOSE_PrintErrors);
    }
    /* Convert to Torr */
    if (pressureMultiplier!=1) {
      long j;
      for (j=0; j<pressureData->nLocations; j++) {
        pressureData->pressure[i][j] *= pressureMultiplier;
      }
    }
  }

  dsMax = -(dsMin = DBL_MAX);
  for (i=1; i<pressureData->nLocations; i++) {
    ds = pressureData->s[i] - pressureData->s[i-1];
    if (ds<=0)
      bombVA("s data is not monotonically increasing in pressure data file %s (%le, %le)", filename, pressureData->s[i-1], pressureData->s[i]);
    if (dsMin>ds)
      dsMin = ds;
    if (dsMax<ds)
      dsMax = ds;
  }
  if (fabs(1-dsMin/dsMax)>1e-3)
      bombVA("s data is not uniformly spaced to within desired 0.1% in pressure data file %s", filename);

}

void computeAverageGasPressures(double sStart, double sEnd, double *pressure, PRESSURE_DATA *pressureData)
{
  double sum;
  long iGas, iLocation, iStart, iEnd;

  /* Could improve this by interpolating the pressure at sStart and sEnd */

  /* Find the indices spanning the desired region */
  iStart = iEnd = -1;
  for (iLocation=0; iLocation<pressureData->nLocations; iLocation++) {
    if (pressureData->s[iLocation]>=sStart && iStart==-1)
      iStart = iLocation;
    if (pressureData->s[iLocation]<=sEnd)
      iEnd = iLocation;
    else
      break;
  }
  if (iStart==-1 || iEnd==-1 || iEnd<=iStart)
    bombVA("Failed to find indices corresponding to pressure region s:[%le, %le] m\n",
                  sStart, sEnd);

  for (iGas=0; iGas<pressureData->nGasses; iGas++) {
    sum = 0;
    for (iLocation=iStart; iLocation<=iEnd; iLocation++) 
      sum += pressureData->pressure[iGas][iLocation];
    pressure[iGas] = sum/(iEnd-iStart+1);
  }
}

#define GAS_H2 0
#define GAS_CH4 1
#define GAS_CO2 2
#define GAS_CO 3
#define GAS_N2 4
#define GAS_H2O 5
#define N_GASSES 6
static char *gasName[N_GASSES] = {"H2", "CH4", "CO2", "CO", "N2", "H2O"};

long identifyGas(GAS_DATA *gasData, char *name)
{
  switch (match_string(name, gasName, N_GASSES, EXACT_MATCH)) {
  case GAS_H2:
    gasData->nConstituents = 1;
    gasData->nAtoms = tmalloc(sizeof(*(gasData->nAtoms))*(gasData->nConstituents));
    gasData->Z = tmalloc(sizeof(*(gasData->Z))*(gasData->nConstituents));
    gasData->A = tmalloc(sizeof(*(gasData->A))*(gasData->nConstituents));
    gasData->nAtoms[0] = 2;
    gasData->Z[0] = 1;
    gasData->A[0] = 1;
    break;
  case GAS_H2O:
    gasData->nConstituents = 2;
    gasData->nAtoms = tmalloc(sizeof(*(gasData->nAtoms))*(gasData->nConstituents));
    gasData->Z = tmalloc(sizeof(*(gasData->Z))*(gasData->nConstituents));
    gasData->A = tmalloc(sizeof(*(gasData->A))*(gasData->nConstituents));
    gasData->nAtoms[0] = 2;
    gasData->Z[0] = 1;
    gasData->A[0] = 1;
    gasData->nAtoms[1] = 1;
    gasData->Z[1] = 8;
    gasData->A[1] = 16;
    break;
  case GAS_CH4:
    gasData->nConstituents = 2;
    gasData->nAtoms = tmalloc(sizeof(*(gasData->nAtoms))*(gasData->nConstituents));
    gasData->Z = tmalloc(sizeof(*(gasData->Z))*(gasData->nConstituents));
    gasData->A = tmalloc(sizeof(*(gasData->A))*(gasData->nConstituents));
    gasData->nAtoms[0] = 1;
    gasData->Z[0] = 6;
    gasData->A[0] = 12;
    gasData->nAtoms[1] = 4;
    gasData->Z[1] = 1;
    gasData->A[1] = 1;
    break;
  case GAS_CO2:
    gasData->nConstituents = 2;
    gasData->nAtoms = tmalloc(sizeof(*(gasData->nAtoms))*(gasData->nConstituents));
    gasData->Z = tmalloc(sizeof(*(gasData->Z))*(gasData->nConstituents));
    gasData->A = tmalloc(sizeof(*(gasData->A))*(gasData->nConstituents));
    gasData->nAtoms[0] = 1;
    gasData->Z[0] = 6;
    gasData->A[0] = 12;
    gasData->nAtoms[1] = 2;
    gasData->Z[1] = 8;
    gasData->A[1] = 16;
    break;
  case GAS_CO:
    gasData->nConstituents = 2;
    gasData->nAtoms = tmalloc(sizeof(*(gasData->nAtoms))*(gasData->nConstituents));
    gasData->Z = tmalloc(sizeof(*(gasData->Z))*(gasData->nConstituents));
    gasData->A = tmalloc(sizeof(*(gasData->A))*(gasData->nConstituents));
    gasData->nAtoms[0] = 1;
    gasData->Z[0] = 6;
    gasData->A[0] = 12;
    gasData->nAtoms[1] = 1;
    gasData->Z[1] = 8;
    gasData->A[1] = 16;
    break;
  case GAS_N2:
    gasData->nConstituents = 1;
    gasData->nAtoms = tmalloc(sizeof(*(gasData->nAtoms))*(gasData->nConstituents));
    gasData->Z = tmalloc(sizeof(*(gasData->Z))*(gasData->nConstituents));
    gasData->A = tmalloc(sizeof(*(gasData->A))*(gasData->nConstituents));
    gasData->nAtoms[0] = 2;
    gasData->Z[0] = 7;
    gasData->A[0] = 14;
    break;
  default:
    return 0;
    break;    
  }
  return 1;
}
