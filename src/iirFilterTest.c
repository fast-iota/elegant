/*************************************************************************\
* Copyright (c) 2015 The University of Chicago, as Operator of Argonne
* National Laboratory.
* Copyright (c) 2015 The Regents of the University of California, as
* Operator of Los Alamos National Laboratory.
* This file is distributed subject to a Software License Agreement found
* in the file LICENSE that is included with this distribution. 
\*************************************************************************/

#include "mdb.h"
#include "track.h"
void bombElegantVA(char *template, ...);

int main(int argc, char **argv) {
  char *input, *output, *filterFile;
  SDDS_DATASET SDDSin, SDDSout;
  IIRFILTER filterBank[4];
  long i, rows, nFilters;
  double *x, *y;

  if (argc < 4)
    bomb("usage: iirFilter <input> <output> <filterFile>", NULL);
  input = argv[1];
  output = argv[2];
  filterFile = argv[3];

  for (i = 0; i < 4; i++)
    memset(&filterBank[i], 0, sizeof(filterBank[i]));

  rows = 0;
  if (!SDDS_InitializeInput(&SDDSin, input) ||
      SDDS_ReadPage(&SDDSin) != 1 ||
      (rows = SDDS_RowCount(&SDDSin)) <= 0)
    SDDS_Bomb("Problem with input file");
  if (!(x = SDDS_GetColumnInDoubles(&SDDSin, "x")) ||
      !SDDS_Terminate(&SDDSin))
    SDDS_Bomb("Problem with input file");

  if (!SDDS_InitializeOutputElegant(&SDDSout, SDDS_BINARY, 1, NULL, NULL, output) ||
      SDDS_DefineColumn(&SDDSout, "y", NULL, NULL, NULL, NULL, SDDS_DOUBLE, 0) < 0 ||
      !SDDS_WriteLayout(&SDDSout) ||
      !SDDS_StartPage(&SDDSout, rows))
    SDDS_Bomb("Problem with output file");

  nFilters = readIIRFilter(filterBank, 4, filterFile);
  y = tmalloc(sizeof(*y) * rows);
  for (i = 0; i < rows; i++)
    y[i] = applyIIRFilter(filterBank, nFilters, x[i]);
  if (!SDDS_SetColumn(&SDDSout, SDDS_SET_BY_NAME, y, rows, "y") || !SDDS_WritePage(&SDDSout) || !SDDS_Terminate(&SDDSout))
    SDDS_Bomb("Problem with output file (2)");

  free(x);
  free(y);

  freeIIRFilterMemory(filterBank, nFilters);

  return 0;
}
