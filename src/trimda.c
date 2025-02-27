/* Copyright 2019 by Michael Borland and Argonne National Laboratory,
 * all rights reserved.
 */

#include "mdb.h"
#include "SDDS.h"
#include "track.h"

static char *USAGE = "trimda <input> <output>";

int main(int argc, char **argv) {
  char *output, *input;
  SDDS_DATASET SDDS_input, SDDS_output;
  long line, lines, full_plane;
  double dtheta;
  double *xLimit, *yLimit;
  double *dxFactor, *dyFactor, area;

  if (argc != 3) {
    fprintf(stderr, "%s\n", USAGE);
    exit(1);
  }
  input = argv[1];
  output = argv[2];

  if (!fexists(input)) {
    fprintf(stderr, "not found: %s\n", input);
    exit(1);
  }
  if (fexists(output)) {
    fprintf(stderr, "in use: %s\n", output);
    exit(1);
  }

  if (!SDDS_InitializeInput(&SDDS_input, input)) {
    fprintf(stderr, "unable to initialize input file\n");
    exit(1);
  }
  if (SDDS_CheckColumn(&SDDS_input, "x", "m", SDDS_ANY_FLOATING_TYPE, stderr) != SDDS_CHECK_OKAY ||
      SDDS_CheckColumn(&SDDS_input, "y", "m", SDDS_ANY_FLOATING_TYPE, stderr) != SDDS_CHECK_OKAY) {
    fprintf(stderr, "missing data or wrong units in input file\n");
    exit(1);
  }
  if (!SDDS_InitializeOutput(&SDDS_output, SDDS_BINARY, 1, NULL, NULL, output) ||
      !SDDS_TransferAllParameterDefinitions(&SDDS_output, &SDDS_input, SDDS_TRANSFER_OVERWRITE) ||
      !SDDS_TransferAllColumnDefinitions(&SDDS_output, &SDDS_input, SDDS_TRANSFER_OVERWRITE) ||
      !SDDS_WriteLayout(&SDDS_output)) {
    SDDS_SetError("Probem setting up output file");
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors | SDDS_EXIT_PrintErrors);
  }

  while (SDDS_ReadPage(&SDDS_input) > 0) {
    lines = SDDS_RowCount(&SDDS_input);
    if (!(xLimit = SDDS_GetColumnInDoubles(&SDDS_input, "x")) ||
        !(yLimit = SDDS_GetColumnInDoubles(&SDDS_input, "y"))) {
      fprintf(stderr, "unable to get x or y data from input file\n");
      exit(1);
    }
    if (!SDDS_StartPage(&SDDS_output, lines)) {
      SDDS_SetError("Probem starting page in output file");
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors | SDDS_EXIT_PrintErrors);
    }
    if (!SDDS_CopyColumns(&SDDS_output, &SDDS_input)) {
      SDDS_SetError("Probem copying column data to output file");
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors | SDDS_EXIT_PrintErrors);
    }
    if (!SDDS_CopyParameters(&SDDS_output, &SDDS_input)) {
      SDDS_SetError("Probem copying parameter data to output file");
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors | SDDS_EXIT_PrintErrors);
    }

    full_plane = 0;
    for (line = 0; line < lines; line++) {
      if (yLimit[line] < 0) {
        full_plane = 1;
        break;
      }
    }

    dxFactor = tmalloc(sizeof(*dxFactor) * lines);
    dyFactor = tmalloc(sizeof(*dyFactor) * lines);
    if (full_plane == 0)
      dtheta = PI / (lines - 1);
    else
      dtheta = PIx2 / (lines - 1);
    for (line = 0; line < lines; line++) {
      dxFactor[line] = sin(-PI / 2 + dtheta * line);
      dyFactor[line] = cos(-PI / 2 + dtheta * line);
      if (fabs(dxFactor[line]) < 1e-6)
        dxFactor[line] = 0;
      if (fabs(dyFactor[line]) < 1e-6)
        dyFactor[line] = 0;
    }

    area = trimApertureSearchResult(lines, xLimit, yLimit, dxFactor, dyFactor, full_plane);

    if (!SDDS_SetColumn(&SDDS_output, SDDS_SET_BY_NAME, xLimit, lines, "xClipped") ||
        !SDDS_SetColumn(&SDDS_output, SDDS_SET_BY_NAME, yLimit, lines, "yClipped") ||
        !SDDS_SetParameters(&SDDS_output, SDDS_SET_BY_NAME | SDDS_PASS_BY_VALUE, "Area", area, NULL)) {
      SDDS_SetError("Probem setting data in output file");
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors | SDDS_EXIT_PrintErrors);
    }
    if (!SDDS_WritePage(&SDDS_output)) {
      SDDS_SetError("Probem writing page to output file");
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors | SDDS_EXIT_PrintErrors);
    }

    free(dxFactor);
    free(dyFactor);
    free(xLimit);
    free(yLimit);
  }
  if (!SDDS_Terminate(&SDDS_input) || !SDDS_Terminate(&SDDS_output)) {
    SDDS_SetError("Probem terminating input and output file");
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors | SDDS_EXIT_PrintErrors);
  }

  exit(0);
}
