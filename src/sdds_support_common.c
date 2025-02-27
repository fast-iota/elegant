/*************************************************************************\
* Copyright (c) 2002 The University of Chicago, as Operator of Argonne
* National Laboratory.
* Copyright (c) 2002 The Regents of the University of California, as
* Operator of Los Alamos National Laboratory.
* This file is distributed subject to a Software License Agreement found
* in the file LICENSE that is included with this distribution. 
\*************************************************************************/

/* file: sdds_support.c
 * purpose: routines for working with SDDS files
 * 
 * Michael Borland, 1993.
 */
#include "mdb.h"
#include "fftpackC.h"
#include "track.h"
#include "SDDS.h"

void setUpFtableBookn(ntuple *book,
                      double *B,
                      double xmin, double xmax, double dx, long nx,
                      double ymin, double ymax, double dy, long ny,
                      double zmin, double zmax, double dz, long nz);

long check_sdds_column(SDDS_TABLE *SDDS_table, char *name, char *units) {
  char *units1;
  if (SDDS_GetColumnIndex(SDDS_table, name) < 0)
    return (0);
  if (SDDS_GetColumnInformation(SDDS_table, "units", &units1, SDDS_GET_BY_NAME, name) != SDDS_STRING) {
    SDDS_SetError("units field of column has wrong data type!");
    SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors | SDDS_VERBOSE_PrintErrors);
  }
  if (!units || SDDS_StringIsBlank(units)) {
    if (!units1)
      return (1);
    if (SDDS_StringIsBlank(units1)) {
      free(units1);
      return (1);
    }
    return (0);
  }
  if (!units1)
    return (0);
  if (strcmp(units, units1) == 0) {
    free(units1);
    return (1);
  }
  free(units1);
  return (0);
}

long check_sdds_parameter(SDDS_TABLE *SDDS_table, char *name, char *units) {
  char *units1;
  if (SDDS_GetParameterIndex(SDDS_table, name) < 0)
    return (0);
  if (SDDS_GetParameterInformation(SDDS_table, "units", &units1, SDDS_GET_BY_NAME, name) != SDDS_STRING) {
    SDDS_SetError("units field of parameter has wrong data type!");
    SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors | SDDS_VERBOSE_PrintErrors);
  }
  if (!units || SDDS_StringIsBlank(units)) {
    if (!units1)
      return (1);
    if (SDDS_StringIsBlank(units1)) {
      free(units1);
      return (1);
    }
    return (0);
  }
  if (!units1)
    return (0);
  if (strcmp(units, units1) == 0) {
    free(units1);
    return (1);
  }
  free(units1);
  return (0);
}

void reorganizeFTABLE(double **F, double **Ftmp, long nx, long ny, long nz);

void readSimpleFtable(FTABLE *ftable) {
  /* Assume input file is same as for BMAPXYZ: (x, y, z, Bx, By, Bz) */
  SDDS_DATASET SDDSin;
  double *x = NULL, *y = NULL, *z = NULL, *Fx = NULL, *Fy = NULL, *Fz = NULL, *Ftmp;
  long nx, ny, nz, points;
  double dx, dy, dz;
  double xmin, ymin, zmin;
  double xmax, ymax, zmax;

  if (!fexists(ftable->inputFile)) {
    printf("file %s not found for FTABLE element\n", ftable->inputFile);
    fflush(stdout);
    exitElegant(1);
  }

  if (!SDDS_InitializeInputFromSearchPath(&SDDSin, ftable->inputFile) ||
      SDDS_ReadPage(&SDDSin) <= 0 ||
      !(x = SDDS_GetColumnInDoubles(&SDDSin, "x")) || !(y = SDDS_GetColumnInDoubles(&SDDSin, "y")) ||
      !(z = SDDS_GetColumnInDoubles(&SDDSin, "z"))) {
    SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors | SDDS_VERBOSE_PrintErrors);
  }
  if (!check_sdds_column(&SDDSin, "x", "m") ||
      !check_sdds_column(&SDDSin, "y", "m") ||
      !check_sdds_column(&SDDSin, "z", "m")) {
    fprintf(stderr, "FTABLE input file must have x, y, and z in m (meters)\n");
    exitElegant(1);
  }

  if (!(Fx = SDDS_GetColumnInDoubles(&SDDSin, "Bx")) || !(Fy = SDDS_GetColumnInDoubles(&SDDSin, "By")) ||
      !(Fz = SDDS_GetColumnInDoubles(&SDDSin, "Bz"))) {
    fprintf(stderr, "FTABLE input file must have (Bx, By, Bz) (in T)\n");
    exitElegant(1);
  }

  if (!check_sdds_column(&SDDSin, "Bx", "T") ||
      !check_sdds_column(&SDDSin, "By", "T") ||
      !check_sdds_column(&SDDSin, "Bz", "T")) {
    fprintf(stderr, "FTABLE input file must have Bx, By, and Bz in T (Tesla)\n");
    exitElegant(1);
  }

  if (!(points = SDDS_CountRowsOfInterest(&SDDSin)) || points < 2) {
    printf("file %s for FTABLE element has insufficient data\n", ftable->inputFile);
    fflush(stdout);
    exitElegant(1);
  }
  SDDS_Terminate(&SDDSin);

  /* It is required that the data is ordered so that x changes fastest.
     * This can be accomplished with sddssort -column=z,incr -column=y,incr -column=x,incr
     * The points are assumed to be equipspaced.
     */
  nx = 1;
  xmin = x[0];
  while (nx < points) {
    if (x[nx - 1] > x[nx])
      break;
    nx++;
  }
  if (nx == points) {
    printf("file %s for FTABLE element doesn't have correct structure or amount of data (x)\n",
           ftable->inputFile);
    printf("Use sddssort -column=z -column=y -column=x to sort the file\n");
    exitElegant(1);
  }
  xmax = x[nx - 1];
  dx = (xmax - xmin) / (nx - 1);

  ny = 1;
  ymin = y[0];
  while (ny < (points / nx)) {
    if (y[(ny - 1) * nx] > y[ny * nx])
      break;
    ny++;
  }
  if (ny == points) {
    printf("file %s for FTABLE element doesn't have correct structure or amount of data (y)\n",
           ftable->inputFile);
    printf("Use sddssort -column=z -column=y -column=x to sort the file\n");
    exitElegant(1);
  }
  ymax = y[(ny - 1) * nx];
  dy = (ymax - ymin) / (ny - 1);
  nz = 0; /* prevent compiler warning */
  if (nx <= 1 || ny <= 1 || (nz = points / (nx * ny)) <= 1) {
    printf("file %s for FTABLE element doesn't have correct structure or amount of data\n",
           ftable->inputFile);
    printf("nx = %ld, ny=%ld, nz=%ld, points=%ld\n", nx, ny, nz, points);
    exitElegant(1);
  }
  zmin = z[0];
  zmax = z[points - 1];
  dz = (zmax - zmin) / (nz - 1);
  printf("FTABLE element from file %s: nx=%ld, ny=%ld, nz=%ld\ndx=%e, dy=%e, dz=%e\nx:[%e, %e], y:[%e, %e], z:[%e, %e]\n",
         ftable->inputFile,
         nx, ny, nz,
         dx, dy, dz,
         xmin, xmax,
         ymin, ymax,
         zmin, zmax);
  free(x);
  free(y);
  free(z);

  Ftmp = tmalloc(sizeof(*Ftmp) * nx * ny * nz);
  reorganizeFTABLE(&Fx, &Ftmp, nx, ny, nz);
  reorganizeFTABLE(&Fy, &Ftmp, nx, ny, nz);
  reorganizeFTABLE(&Fz, &Ftmp, nx, ny, nz);
  free(Ftmp);

  ftable->Bx = (ntuple *)malloc(sizeof(*(ftable->Bx)));
  ftable->By = (ntuple *)malloc(sizeof(*(ftable->By)));
  ftable->Bz = (ntuple *)malloc(sizeof(*(ftable->Bz)));

  setUpFtableBookn(ftable->Bx, Fx, xmin, xmax, dx, nx, ymin, ymax, dy, ny, zmin, zmax, dz, nz);
  setUpFtableBookn(ftable->By, Fy, xmin, xmax, dx, nx, ymin, ymax, dy, ny, zmin, zmax, dz, nz);
  setUpFtableBookn(ftable->Bz, Fz, xmin, xmax, dx, nx, ymin, ymax, dy, ny, zmin, zmax, dz, nz);
}

void setUpFtableBookn(ntuple *book,
                      double *B,
                      double xmin, double xmax, double dx, long nx,
                      double ymin, double ymax, double dy, long ny,
                      double zmin, double zmax, double dz, long nz) {
  long i;

  book->nD = 3;

  book->vname = calloc(sizeof(*book->vname), book->nD);
  book->units = calloc(sizeof(*book->units), book->nD);
  book->xmin = calloc(sizeof(*book->xmin), book->nD);
  book->xmax = calloc(sizeof(*book->xmax), book->nD);
  book->dx = calloc(sizeof(*book->dx), book->nD);
  book->xbins = calloc(sizeof(*book->xbins), book->nD);

  book->xmin[0] = xmin;
  book->xmax[0] = xmax;
  book->dx[0] = dx;
  book->xbins[0] = nx;

  book->xmin[1] = ymin;
  book->xmax[1] = ymax;
  book->dx[1] = dy;
  book->xbins[1] = ny;

  book->xmin[2] = zmin;
  book->xmax[2] = zmax;
  book->dx[2] = dz;
  book->xbins[2] = nz;

  book->length = nx * ny * nz;
  book->total = 0;
  book->value = B;
  for (book->total = 0, i = 0; i < book->length; i++)
    book->total += book->value[i];
}

/* The simple FTABLE input file is ordered differently (to conform with BMAPXYZ), so we have
 * to re-order it.
 */

void reorganizeFTABLE(double **F, double **Ftmp, long nx, long ny, long nz) {
  long ix, iy, iz;
  long j, k;
  double *Fsave;
  for (ix = 0; ix < nx; ix++) {
    for (iy = 0; iy < ny; iy++) {
      for (iz = 0; iz < nz; iz++) {
        j = ix + iy * nx + iz * nx * ny;
        k = iz + iy * nz + ix * nz * ny;
        (*Ftmp)[k] = (*F)[j];
      }
    }
  }
  Fsave = *F;
  *F = *Ftmp;
  *Ftmp = Fsave;
}

int32_t SDDS_InitializeOutputElegant(SDDS_DATASET *SDDS_dataset, int32_t data_mode,
                                     int32_t lines_per_row, const char *description,
                                     const char *contents, const char *filename) {
  char *s;
  if (allowOverwrite || !fexists(filename))
    return SDDS_InitializeOutput(SDDS_dataset, data_mode, lines_per_row, description,
                                 contents, filename);
  s = tmalloc(sizeof(*filename) * (strlen(filename) + 1000));
  sprintf(s, "Unable to write to %s --- file exists", filename);
  SDDS_SetError(s);
  return 0;
}

int32_t SDDS_Parallel_InitializeOutputElegant(SDDS_DATASET *SDDS_dataset, const char *description, const char *contents, const char *filename) {
  char *s;
  if (allowOverwrite || !fexists(filename))
    return SDDS_Parallel_InitializeOutput(SDDS_dataset, description, contents, filename);
  s = tmalloc(sizeof(*filename) * (strlen(filename) + 1000));
  sprintf(s, "Unable to write to %s --- file exists", filename);
  SDDS_SetError(s);
  return 0;
}
