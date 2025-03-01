#include "mdb.h"
#include "SDDS.h"
#include "scan.h"
#include "fftpackC.h"

#if defined(linux) || (defined(_WIN32) && !defined(_MINGW))
#  include <omp.h>
#else
#  define NOTHREADS 1
#endif

typedef struct {
  double *x, dx, *z, dz;
  double **By;
  long nx, nz;
  long nxMax;
} MIDPLANE_FIELDS;

typedef struct {
  double dx, dz;
  long nx, nz;
  double **DzBy, **DxBy;
  double **Dx2By, **Dz2By;
  double **Dz3By, **Dx3By;
  double **DzDx2By, **Dz2DxBy;
  double **Dx4By, **Dz2Dx2By, **Dz4By;
  double **Dz5By, **Dz3Dx2By, **DzDx4By;
  double **Dx5By, **Dz2Dx3By, **Dz4DxBy;
} MIDPLANE_EXPANSION;

int ReadInputFile(MIDPLANE_FIELDS *mpFields, char *xParameterName, char *zColumnName, char *ByColumnName, char *input);
void ComputeDerivatives(MIDPLANE_EXPANSION *expansion, MIDPLANE_FIELDS *mpFields,
                        long xDerivOrder, long zDerivOrder);
int WriteExpansion(char *output, MIDPLANE_EXPANSION *expansion);
int EvaluateExpansionAndOutput(char *output, long ny, double yMax, MIDPLANE_EXPANSION *expansion, MIDPLANE_FIELDS *fields);
void SmoothData(MIDPLANE_FIELDS *fields, long smoothPoints, long fitTerms, double fraction, short versus_x, char *outputFile);

#define SET_YORDER 0
#define SET_EVALUATE 1
#define SET_INPUT 2
#define SET_OUTPUT 3
#define SET_DERIVATIVE_ORDER 4
#define SET_X_SMOOTHING 5
#define SET_Z_SMOOTHING 6
#define SET_SMOOTHING_PASSES 7
#define N_OPTIONS 8

char *option[N_OPTIONS] = {
  "yorder", "evaluate", "input", "output", "derivativeorder", "xsmoothing", "zsmoothing", "smoothingpasses"
};

#define USAGE "offMidplaneExpansion -input=<filename>[,xName=<parameterName>][,zName=<columnName>][,ByName=<columnName>]\n\
              [-yOrder=<integer>] [-derivativeOrder=z={2|4|6},x={2|4|6}]\n\
              [-xSmoothing=[points=<integer>,][terms=<integer>][,output=<filename>][,fraction=<value>]]\n\
              [-zSmoothing=[points=<integer>,][terms=<integer>][,output=<filename>][,fraction=<value>]\n\
              [-smoothingPasses=<number>]\n\
              [-evaluate=<filename>,ny=<integer>,yhalfspan=<meters>]\n\
-input       Provides the name of an SDDS file containing By on a grid of (x, z) points.\n\
             The file must have one page for each value of x, with the value of x\n\
             given in the named parameter. The z values must be identical on all pages.\n\
-yOrder      Order of the expansion in the vertical coordinate. Defaults to 3.\n\
-derivativeOrder\n\
             Allows setting the order of the finite-difference equations used for taking\n\
             derivatives versus z and x. E.g., order 2 means that the error is O(h^2).\n\
             Higher order will amplify noise but provide better reproduction of sharp\n\
             features. Defaults to 4.\n\
-xSmoothing\n\
-zSmoothing\n\
             If given, then for the indicated plane, data are smoothed using Savitzky-Golay\n\
             filters prior to taking derivatives. The filters use least-squares fits over the\n\
             given number of data points. By default, 5 points are used with a quadratic (3-term) fit.\n\
             The number of points must be odd. The fraction parameter specifies what fraction of\n\
             the smoothed result should be blended with the original data; it defaults to 1, which means\n\
             the original data is replaced; a value of 0 means the original data is used unaltered.\n\
-evaluate    Evaluates the expansion on a 3D grid. This file can be used with the BRAT element in\n\
             ELEGANT.\n\n\
Program by Michael Borland, 2023."

int main(int argc, char **argv)
{
  SCANNED_ARG *scanned;
  long i_arg, yOrder = 3;
  int32_t xDerivOrder = 4, zDerivOrder = 4;
  char *input = NULL, *output = NULL, *evaluationOutput = NULL;
  char *xParameter = "x", *zColumn = "z", *ByColumn = "By";
  unsigned long dummyFlags;
  int32_t evaluationPoints = 0;
  double evaluationYMax = 0;
  MIDPLANE_FIELDS midplaneFields;
  MIDPLANE_EXPANSION midplaneExpansion;
  int32_t xFitPoints, zFitPoints, xFitTerms, zFitTerms;
  char *xSmoothedOutput, *zSmoothedOutput;
  double xSmoothingFraction, zSmoothingFraction;
  long smoothingPasses = 1;

  argc = scanargs(&scanned, argc, argv);
  if (argc < 2 || argc > (2 + N_OPTIONS)) {
    fprintf(stderr, "%s\n", USAGE);
    return (1);
  }
  
  xFitPoints = zFitPoints = xFitTerms = zFitTerms = 0;
  xSmoothedOutput = zSmoothedOutput = NULL;
  xSmoothingFraction = zSmoothingFraction = 1;

  for (i_arg = 1; i_arg < argc; i_arg++) {
    if (scanned[i_arg].arg_type == OPTION) {
      /* process options here */
      switch (match_string(scanned[i_arg].list[0], option, N_OPTIONS, 0)) {
      case SET_INPUT:
        if (scanned[i_arg].n_items<2 || scanned[i_arg].n_items>5) {
          fprintf(stderr, "invalid -input syntax\n%s\n", USAGE);
          return (1);
        }
        input = scanned[i_arg].list[1];
        scanned[i_arg].n_items -= 2;
        if (!scanItemList(&dummyFlags, scanned[i_arg].list+2, &scanned[i_arg].n_items, 0,
                          "xname", SDDS_STRING, &xParameter, 1, 0, 
                          "zname", SDDS_STRING, &zColumn, 1, 0, 
                          "byname", SDDS_STRING, &ByColumn, 1, 0, 
                          NULL)) {
          fprintf(stderr, "invalid -input syntax\n%s\n", USAGE);
          return (1);
        }
        break;
      case SET_OUTPUT:
        if (scanned[i_arg].n_items!=2 || scanned[i_arg].n_items>3) {
          fprintf(stderr, "invalid -output syntax\n%s\n", USAGE);
          return (1);
        }
        output = scanned[i_arg].list[1];
        break;
      case SET_EVALUATE:
        if (scanned[i_arg].n_items < 4) {
          fprintf(stderr, "invalid -evaluate syntax\n%s\n", USAGE);
          return (1);
        }
        evaluationOutput = scanned[i_arg].list[1];
        scanned[i_arg].n_items -= 2;
        evaluationPoints = -1;
        evaluationYMax = 0;
        if (scanned[i_arg].n_items>0 &&
            !scanItemList(&dummyFlags, scanned[i_arg].list+2, &scanned[i_arg].n_items, 0,
                          "yhalfspan", SDDS_DOUBLE, &evaluationYMax, 1, 0,
                          "ny", SDDS_LONG, &evaluationPoints, 1, 0,
                          NULL)) {
          fprintf(stderr, "invalid -evaluation syntax\n%s\n", USAGE);
          exit(1);
        }
        if (evaluationPoints<2)
          SDDS_Bomb("ny must be 2 or greater");
        if (evaluationYMax<=0)
          SDDS_Bomb("yhalfspan must be positive");
        break;
      case SET_YORDER:
        if (scanned[i_arg].n_items != 2 || sscanf(scanned[i_arg].list[1], "%ld", &yOrder) != 1 || yOrder < 1)
          SDDS_Bomb("invalid -yOrder syntax: give an value greater than 0");
        break;
      case SET_DERIVATIVE_ORDER:
        if (scanned[i_arg].n_items <2 || scanned[i_arg].n_items>4) {
          fprintf(stderr, "Invalid -derivativeOrder syntax\n%s\n", USAGE);
          exit(1);
        }
        scanned[i_arg].n_items -= 1;
        if (!scanItemList(&dummyFlags, scanned[i_arg].list+1, &scanned[i_arg].n_items, 0,
                          "x", SDDS_LONG, &xDerivOrder, 1, 0,
                          "z", SDDS_LONG, &zDerivOrder, 1, 0,
                          NULL) || 
            (xDerivOrder!=2 && xDerivOrder!=4 && xDerivOrder!=6) ||
            (zDerivOrder!=2 && zDerivOrder!=4 && zDerivOrder!=6) ) {
          fprintf(stderr, "invalid -derivativeOrder syntax\n%s\n", USAGE);
          exit(1);
        }
        break;
      case SET_X_SMOOTHING:
        xFitPoints = 5;
        xFitTerms = 3;
        xSmoothingFraction = 1;
        xSmoothedOutput = NULL;
        if ((scanned[i_arg].n_items -= 1)>0) {
          if (!scanItemList(&dummyFlags, scanned[i_arg].list+1, &scanned[i_arg].n_items, 0,
                            "points", SDDS_LONG, &xFitPoints, 1, 0,
                            "terms", SDDS_LONG, &xFitTerms, 1, 0,
                            "output", SDDS_STRING, &xSmoothedOutput, 1, 0,
                            "fraction", SDDS_DOUBLE, &xSmoothingFraction, 1, 0,
                            NULL) || 
              (xFitPoints!=0 && (xFitPoints<3 || xFitPoints%2!=1 || xFitTerms>xFitPoints)) || xSmoothingFraction<=0 ||
              xSmoothingFraction>1 ) {
            fprintf(stderr, "invalid -xSmoothing syntax\n%s\n", USAGE);
            exit(1);
          }
        }
        break;
      case SET_Z_SMOOTHING:
        zFitPoints = 5;
        zFitTerms = 3;
        zSmoothingFraction = 1;
        zSmoothedOutput = NULL;
        if ((scanned[i_arg].n_items -= 1)>0) {
          if (!scanItemList(&dummyFlags, scanned[i_arg].list+1, &scanned[i_arg].n_items, 0,
                            "points", SDDS_LONG, &zFitPoints, 1, 0,
                            "terms", SDDS_LONG, &zFitTerms, 1, 0,
                            "output", SDDS_STRING, &zSmoothedOutput, 1, 0,
                            "fraction", SDDS_DOUBLE, &zSmoothingFraction, 1, 0,
                            NULL) || 
              (zFitPoints!=0 && (zFitPoints<3 || zFitPoints%2!=1 || zFitTerms>zFitPoints)) || zSmoothingFraction<=0 ||
              zSmoothingFraction>1 ) {
            fprintf(stderr, "invalid -zSmoothing syntax\n%s\n", USAGE);
            exit(1);
          }
        }
        break;
      case SET_SMOOTHING_PASSES:
        if (scanned[i_arg].n_items != 2 || sscanf(scanned[i_arg].list[1], "%ld", &smoothingPasses) != 1 || smoothingPasses < 1)
          SDDS_Bomb("invalid -smoothingPasses syntax: give an value greater than 0");
        break;
      default:
        fprintf(stderr, "unknown option given\n%s\n", USAGE);
        return (1);
        break;
      }
    } else
      SDDS_Bomb("unrecognized argument");
  }
  
  if (input==NULL) {
    fprintf(stderr, "%s\n", USAGE);
    return (1);
  }
  if (output!=NULL) {
    fprintf(stderr, "Error: output generation not implemented. Use -evaluate option instead.\n");
    return 1;
  }
  if (!ReadInputFile(&midplaneFields, xParameter, zColumn, ByColumn, input))
    SDDS_Bomb("unable to read input file");

  while (smoothingPasses--) {
    if (zFitPoints)
      SmoothData(&midplaneFields, zFitPoints, zFitTerms, zSmoothingFraction, 0, zSmoothedOutput);
    if (xFitPoints) 
      SmoothData(&midplaneFields, xFitPoints, xFitTerms, xSmoothingFraction, 1, xSmoothedOutput);
  }

  ComputeDerivatives(&midplaneExpansion, &midplaneFields, xDerivOrder, zDerivOrder);

  if (!WriteExpansion(output, &midplaneExpansion))
    SDDS_Bomb("unable to write output file");

  if (evaluationOutput)
    if (!EvaluateExpansionAndOutput(evaluationOutput, evaluationPoints, evaluationYMax, &midplaneExpansion, &midplaneFields)) 
      SDDS_Bomb("Problem evaluating expansion");

  return (0);
}

int EvaluateExpansionAndOutput(char *output, long ny, double yMax, MIDPLANE_EXPANSION *xp, MIDPLANE_FIELDS *fields)
{
  long ix, iy, iz, nx, nz;
  SDDS_DATASET SDDSout;
  double y, dy;
  
  if (!SDDS_InitializeOutput(&SDDSout, SDDS_BINARY, 1, NULL, NULL, output) ||
      !SDDS_DefineSimpleParameter(&SDDSout, "y", "m", SDDS_DOUBLE) || 
      !SDDS_DefineSimpleParameter(&SDDSout, "dx", "m", SDDS_DOUBLE) || 
      !SDDS_DefineSimpleParameter(&SDDSout, "dz", "m", SDDS_DOUBLE) || 
      !SDDS_DefineSimpleColumn(&SDDSout, "x", "m", SDDS_DOUBLE) || 
      !SDDS_DefineSimpleColumn(&SDDSout, "z", "m", SDDS_DOUBLE) || 
      !SDDS_DefineSimpleColumn(&SDDSout, "Bx0", "T", SDDS_DOUBLE) || 
      !SDDS_DefineSimpleColumn(&SDDSout, "Bx1", "T", SDDS_DOUBLE) || 
      !SDDS_DefineSimpleColumn(&SDDSout, "Bx2", "T", SDDS_DOUBLE) || 
      !SDDS_DefineSimpleColumn(&SDDSout, "By0", "T", SDDS_DOUBLE) || 
      !SDDS_DefineSimpleColumn(&SDDSout, "By1", "T", SDDS_DOUBLE) || 
      !SDDS_DefineSimpleColumn(&SDDSout, "By2", "T", SDDS_DOUBLE) || 
      !SDDS_DefineSimpleColumn(&SDDSout, "Bz0", "T", SDDS_DOUBLE) ||
      !SDDS_DefineSimpleColumn(&SDDSout, "Bz1", "T", SDDS_DOUBLE) ||
      !SDDS_DefineSimpleColumn(&SDDSout, "Bz2", "T", SDDS_DOUBLE) ||
      !SDDS_WriteLayout(&SDDSout))
    return 0;

  dy = 2*yMax/(ny-1);
  nz = xp->nz;
  nx = xp->nx;

  for (iy=0; iy<ny; iy++) {
    y = -yMax + iy*dy;
    if (!SDDS_StartPage(&SDDSout, nx*nz) || 
        SDDS_SetParameters(&SDDSout, SDDS_SET_BY_NAME|SDDS_PASS_BY_VALUE, 
                           "y", y, "dx", xp->dx, "dz", xp->dz, NULL)!=1)
      return 0;

    /* By */
    for (ix=0; ix<nx; ix++) {
      for (iz=0; iz<nz; iz++)
        if (!SDDS_SetRowValues(&SDDSout, SDDS_SET_BY_NAME|SDDS_PASS_BY_VALUE,
                               ix*nz+iz, 
                               "x", fields->x[ix], 
                               "z", fields->z[iz], 
                               "By0", fields->By[ix][iz],
                               "By1", fields->By[ix][iz] - sqr(y)/2*(xp->Dz2By[ix][iz] + xp->Dx2By[ix][iz]),
                               "By2", fields->By[ix][iz] - sqr(y)/2*(xp->Dz2By[ix][iz] + xp->Dx2By[ix][iz]) +
                               ipow4(y)/24*(xp->Dx4By[ix][iz] + 2*xp->Dz2Dx2By[ix][iz] + xp->Dz4By[ix][iz]),
                               NULL))
          return 0;
    }

    /* Bx */
    for (ix=0; ix<nx; ix++) {
      for (iz=0; iz<nz; iz++)
        if (!SDDS_SetRowValues(&SDDSout, SDDS_SET_BY_NAME|SDDS_PASS_BY_VALUE,
                               ix*nz+iz, 
                               "Bx0", y*xp->DxBy[ix][iz],
                               "Bx1", y*xp->DxBy[ix][iz] - ipow3(y)/6*(xp->Dz2DxBy[ix][iz] + xp->Dx3By[ix][iz]),
                               "Bx2", y*xp->DxBy[ix][iz] - ipow3(y)/6*(xp->Dz2DxBy[ix][iz] + xp->Dx3By[ix][iz]),
                               + ipow5(y)/120*(xp->Dx5By[ix][iz] + 2*xp->Dz2Dx3By[ix][iz] + xp->Dz4DxBy[ix][iz]),
                               NULL))
          return 0;
    }

    /* Bz */
    for (ix=0; ix<nx; ix++) {
      for (iz=0; iz<nz; iz++) {
        if (!SDDS_SetRowValues(&SDDSout, SDDS_SET_BY_NAME|SDDS_PASS_BY_VALUE,
                               ix*nz+iz, 
                               "Bz0", y*xp->DzBy[ix][iz],
                               "Bz1", y*xp->DzBy[ix][iz] - ipow3(y)/6*(xp->Dz3By[ix][iz] + xp->DzDx2By[ix][iz]),
                               "Bz2", y*xp->DzBy[ix][iz] - ipow3(y)/6*(xp->Dz3By[ix][iz] + xp->DzDx2By[ix][iz]) +
                               ipow5(y)/120*(xp->Dz5By[ix][iz] + 2*xp->Dz3Dx2By[ix][iz] + xp->DzDx4By[ix][iz]),
                               NULL))
          return 0;    
      }
    }

    if (!SDDS_WritePage(&SDDSout))
      return 0;
  }
  return 1;
}

int WriteExpansion(char *output, MIDPLANE_EXPANSION *xp)
{
  return 1;
}

int ReadInputFile(MIDPLANE_FIELDS *mpFields, char *xParameterName, char *zColumnName, char *ByColumnName, char *input)
{
  SDDS_DATASET SDDSin;
  short first;
  long nz, iz, ix;
  double *z, *By, x;
  double dx, dxAve;

  if (!SDDS_InitializeInput(&SDDSin, input)) {
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
    return(0);
  }
  if (SDDS_CheckColumn(&SDDSin, zColumnName, "m", SDDS_ANY_FLOATING_TYPE, stderr) != SDDS_CHECK_OKAY) {
    SDDS_SetError("Named 'z' column nonexistent, wrong type, or wrong units.");
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
    return 0;
  }
  if (SDDS_CheckColumn(&SDDSin, ByColumnName, "T", SDDS_ANY_FLOATING_TYPE, stderr) != SDDS_CHECK_OKAY) {
    SDDS_SetError("Named 'By' column nonexistent, wrong type, or wrong units.");
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
    return 0;
  }
  if (SDDS_CheckParameter(&SDDSin, xParameterName, "m", SDDS_ANY_FLOATING_TYPE, stderr) != SDDS_CHECK_OKAY) {
    SDDS_SetError("Named 'x' parameter nonexistent, wrong type, or wrong units.");
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
    return 0;
  }

  mpFields->nx = mpFields->nz = 0;
  mpFields->x = malloc(sizeof(*(mpFields->x))*(mpFields->nxMax=10));
  mpFields->By = malloc(sizeof(*(mpFields->By))*mpFields->nxMax);

  first = 1;
  while (SDDS_ReadPage(&SDDSin)>0) {
    if ((nz = SDDS_RowCount(&SDDSin))<2) {
      SDDS_SetError("Too few z points");
      return 0;
    }
    if (!first && nz!=mpFields->nz) {
      SDDS_SetError("inconsistent number of rows input file");
      return 0;
    }
    if (!(z=SDDS_GetColumnInDoubles(&SDDSin, zColumnName)) ||
        !(By=SDDS_GetColumnInDoubles(&SDDSin, ByColumnName)) ||
        !SDDS_GetParameterAsDouble(&SDDSin, xParameterName, &x)) {
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
      return(0);
    }
    if (!first) {
      /* check consistency of z values */
      for (iz=0; iz<nz; iz++)
        if (z[iz]!=mpFields->z[iz]) {
          SDDS_SetError("inconsistent z values in input file");
          return 0;
        }
      if (mpFields->x[mpFields->nx-1]>=x) {
        SDDS_SetError("pages are not in order of increasing x");
        return 0;
      }
    } else
      mpFields->nz = nz;
    if (mpFields->nx>=mpFields->nxMax) {
      mpFields->x = SDDS_Realloc(mpFields->x, sizeof(*(mpFields->x))*(mpFields->nxMax+=10));
      mpFields->By = SDDS_Realloc(mpFields->By, sizeof(*(mpFields->By))*mpFields->nxMax);
    }
    mpFields->x[mpFields->nx] = x;
    mpFields->By[mpFields->nx] = By;
    mpFields->nx += 1;
    if (first) {
      double dz, dzAve;
      mpFields->z = z;
      first = 0;
      mpFields->dz = dzAve = (z[nz-1]-z[0])/(nz-1);
      for (iz=1; iz<nz; iz++) {
        dz = z[iz]-z[iz-1];
        if (dz<=0) {
          SDDS_SetError("z values are not monotonically increasing");
          return 0;
        }
        if (fabs(dz-dzAve)/dzAve > 1e-6) {
          SDDS_SetError("z values are not equispaced (deviation greater than 1 ppm)");
          return 0;
        }
      }
    }
  }
  if (mpFields->nx<3) {
    SDDS_SetError("Must have data on at least 3 x planes.");
    return 0;
  }
  if (mpFields->nx%2==0) {
    SDDS_SetError("Number of x planes must be odd.");
    return 0;
  }
  if ((dxAve = (mpFields->x[mpFields->nx-1] - mpFields->x[0])/(mpFields->nx-1))<=0) {
    SDDS_SetError("Average x spacing is not positive");
    return 0;
  }
  for (ix=1; ix<mpFields->nx; ix++) {
    dx = mpFields->x[ix] - mpFields->x[ix-1];
    if (fabs(dx-dxAve)/dxAve > 1e-6) {
      SDDS_SetError("x values are not equispaced (deviation greater than 1 ppm)");
      return 0;
    }
  }
  mpFields->dx = dxAve;
  return 1;
}

/* This is used in an attempt to reduce noise. I don't think it helps. */
double Kahan(long length, double a[], double *error) {

  double sum = 0.0, C = 0.0, Y, T;
  long i;

  /* Kahan's summation formula */
  for (i = 0; i < length; i++) {
    Y = a[i] - C;
    T = sum + Y;
    C = (T - sum) - Y;
    sum = T;
  }

  *error = -C;
  return sum;
}

void takeDerivative(double *y, int64_t rows, double delta, double *deriv, long differenceOrder)
{
  int64_t i;
  double term[7];
  double error;

  if (differenceOrder==2) {
    /* three-point formula */
    /* See https://web.media.mit.edu/~crtaylor/calculator.html */
    for (i=0; i<rows; i++) {
      error = 0;
      if (i==0) {
        term[0] = -3*y[0];
        term[1] = 4*y[1];
        term[2] = -y[2];
      } else if (i==(rows-1)) {
        term[0] = 3*y[i];
        term[1] = -4*y[i-1];
        term[2] = y[i-2];
      } else {
        term[0] = y[i+1];
        term[1] = -y[i-1];
        term[2] = 0;
      }
      deriv[i] = Kahan(3, term, &error)/(2*delta);
    }
  } else if (differenceOrder==4) {
    /* five-point formula --- amplifies noise compared to three-point formula */
    /* See https://web.media.mit.edu/~crtaylor/calculator.html */
    for (i=0; i<rows; i++) {
      error = 0;
      if (i==0) {
        term[0] = -25*y[i];
        term[1] = 48*y[i+1];
        term[2] = -36*y[i+2];
        term[3] = 16*y[i+3];
        term[4] = -3*y[i+4];
      } else if (i==1) {
        term[0] = -3*y[i-1];
        term[1] = -10 *y[i];
        term[2] = 18*y[i+1];
        term[3] = -6*y[i+2];
        term[4] = y[i+3];
      } else if (i==(rows-2)) {
        term[0] = -1*y[i-3];
        term[1] = 6*y[i-2];
        term[2] = -18*y[i-1];
        term[3] = 10*y[i];
        term[4] = 3*y[i+1];
      } else if (i==(rows-1)) {
        term[0] = 3*y[i-4];
        term[1] = -16*y[i-3];
        term[2] = 36*y[i-2];
        term[3] = -48*y[i-1];
        term[4] = 25*y[i];
      } else {
        term[0] = y[i-2];
        term[1] = -8*y[i-1]; 
        term[2] = 8*y[i+1];
        term[3] = -y[i+2];
        term[4] = 0;
      }
      deriv[i] = Kahan(5, term, &error)/(12*delta);
    }
  } else if (differenceOrder==6) {
    /* seven-point formula --- amplifies noise compared to five-point formula */
    /* See https://web.media.mit.edu/~crtaylor/calculator.html */
    for (i=0; i<rows; i++) {
      error = 0;
      if (i==0) {
        term[0] = -147*y[i];
        term[1] = 360*y[i+1];
        term[2] = -450*y[i+2];
        term[3] = 400*y[i+3];
        term[4] = -225*y[i+4];
        term[5] = 72*y[i+5];
        term[6] = -10*y[i+6];
      } else if (i==1) {
        term[0] = -10*y[i-1];
        term[1] = -77*y[i];
        term[2] = 150*y[i+1];
        term[3] = -100*y[i+2];
        term[4] = 50*y[i+3];
        term[5] = -15*y[i+4];
        term[6] = 2*y[i+5];
      } else if (i==2) {
        term[0] = 2*y[i-2];
        term[1] = -24*y[i-1];
        term[2] = -35*y[i];
        term[3] = 80*y[i+1];
        term[4] = -30*y[i+2];
        term[5] = 8*y[i+3];
        term[6] = -y[i+4];
      } else if (i==(rows-3)) {
        term[0] = y[i-4];
        term[1] = -8*y[i-3];
        term[2] = 30*y[i-2];
        term[3] = -80*y[i-1];
        term[4] = 35*y[i];
        term[5] = 24*y[i+1];
        term[6] = -2*y[i+2];
      } else if (i==(rows-2)) {
        term[0] = -2*y[i-5];
        term[1] = 15*y[i-4];
        term[2] = -50*y[i-3];
        term[3] = 100*y[i-2];
        term[4] = -150*y[i-1];
        term[5] = 77*y[i];
        term[6] = 10*y[i+1];
      } else if (i==(rows-1)) {
        term[0] = 10*y[i-6];
        term[1] = -72*y[i-5];
        term[2] = 225*y[i-4];
        term[3] = -400*y[i-3];
        term[4] = 450*y[i-2];
        term[5] = -360*y[i-1];
        term[6] = 147*y[i];
      } else {
        term[0] = -1*y[i-3];
        term[1] = 9*y[i-2];
        term[2] = -45*y[i-1];
        term[3] = 0;
        term[4] = 45*y[i+1];
        term[5] = -9*y[i+2];
        term[6] = y[i+3];
      }
      deriv[i] = Kahan(7, term, &error)/(60*delta);
    }
  } else
    SDDS_Bomb("unknown derivative order");
}

void takeMultipleDerivative(double *y, int64_t rows, double delta, double *deriv, long derivOrder, long differenceOrder) 
{
  double *buffer;

  /* copy the data */
  buffer = malloc(sizeof(*buffer)*rows);
  memcpy(buffer, y, sizeof(*buffer)*rows);

  /* take derivatives */
  while (derivOrder--) {
    takeDerivative(buffer, rows, delta, deriv, differenceOrder);
    memcpy(buffer, deriv, sizeof(*buffer)*rows);
  }
  free(buffer);
}

void ComputeDerivatives(MIDPLANE_EXPANSION *xp, MIDPLANE_FIELDS *fields, long xDerivOrder, long zDerivOrder)
{
  long ix, iz;
  double *input, *deriv;

  xp->dx = fields->dx;
  xp->dz = fields->dz;
  xp->nx = fields->nx;
  xp->nz = fields->nz;

  /* dBy/dz, dBy/dz^2, dBy/dz^3, dBy/dz^4, dBy/dz^5 for each x */
  xp->DzBy  = (double**)czarray_2d(sizeof(double), fields->nx, fields->nz);
  xp->Dz2By = (double**)czarray_2d(sizeof(double), fields->nx, fields->nz);
  xp->Dz3By = (double**)czarray_2d(sizeof(double), fields->nx, fields->nz);
  xp->Dz4By = (double**)czarray_2d(sizeof(double), fields->nx, fields->nz);
  xp->Dz5By = (double**)czarray_2d(sizeof(double), fields->nx, fields->nz);
  for (ix=0; ix<fields->nx; ix++) {
    takeMultipleDerivative(fields->By[ix], fields->nz, fields->dz, xp->DzBy[ix], 1, zDerivOrder);
    takeMultipleDerivative(fields->By[ix], fields->nz, fields->dz, xp->Dz2By[ix], 2, zDerivOrder);
    takeMultipleDerivative(fields->By[ix], fields->nz, fields->dz, xp->Dz3By[ix], 3, zDerivOrder);
    takeMultipleDerivative(fields->By[ix], fields->nz, fields->dz, xp->Dz4By[ix], 4, zDerivOrder);
    takeMultipleDerivative(fields->By[ix], fields->nz, fields->dz, xp->Dz5By[ix], 5, zDerivOrder);
  } 

  /* dBy/dx, dBy/dx^2, dBy/dx^3, ... for each z */
  xp->DxBy  = (double**)czarray_2d(sizeof(double), fields->nx, fields->nz);
  xp->Dx2By = (double**)czarray_2d(sizeof(double), fields->nx, fields->nz);
  xp->Dx3By = (double**)czarray_2d(sizeof(double), fields->nx, fields->nz);
  xp->Dx4By = (double**)czarray_2d(sizeof(double), fields->nx, fields->nz);
  xp->Dx5By = (double**)czarray_2d(sizeof(double), fields->nx, fields->nz);
  input = tmalloc(sizeof(*input)*fields->nx);
  deriv = tmalloc(sizeof(*deriv)*fields->nx);
  for (iz=0; iz<fields->nz; iz++) {
    for (ix=0; ix<fields->nx; ix++)
      input[ix] = fields->By[ix][iz];

    takeMultipleDerivative(input, fields->nx, fields->dx, deriv, 1, xDerivOrder);
    for (ix=0; ix<fields->nx; ix++)
      xp->DxBy[ix][iz] = deriv[ix];

    takeMultipleDerivative(input, fields->nx, fields->dx, deriv, 2, xDerivOrder);
    for (ix=0; ix<fields->nx; ix++)
      xp->Dx2By[ix][iz] = deriv[ix];

    takeMultipleDerivative(input, fields->nx, fields->dx, deriv, 3, xDerivOrder);
    for (ix=0; ix<fields->nx; ix++)
      xp->Dx3By[ix][iz] = deriv[ix];

    takeMultipleDerivative(input, fields->nx, fields->dx, deriv, 4, xDerivOrder);
    for (ix=0; ix<fields->nx; ix++)
      xp->Dx4By[ix][iz] = deriv[ix];

    takeMultipleDerivative(input, fields->nx, fields->dx, deriv, 5, xDerivOrder);
    for (ix=0; ix<fields->nx; ix++)
      xp->Dx5By[ix][iz] = deriv[ix];
  }
  free(input);
  free(deriv);

  /* Dz (Dx^2 By) */
  xp->DzDx2By = (double**)czarray_2d(sizeof(double), fields->nx, fields->nz);
  for (ix=0; ix<fields->nx; ix++) 
    takeMultipleDerivative(xp->Dx2By[ix], fields->nz, fields->dz, xp->DzDx2By[ix], 1, zDerivOrder);

  /* Dz^2 (Dx By) */
  xp->Dz2DxBy = (double**)czarray_2d(sizeof(double), fields->nx, fields->nz);
  for (ix=0; ix<fields->nx; ix++) 
    takeMultipleDerivative(xp->DxBy[ix], fields->nz, fields->dz, xp->Dz2DxBy[ix], 2, zDerivOrder);

  /* Dz^2 (Dx^2 By) */
  xp->Dz2Dx2By = (double**)czarray_2d(sizeof(double), fields->nx, fields->nz);
  for (ix=0; ix<fields->nx; ix++)
    takeMultipleDerivative(xp->Dx2By[ix], fields->nz, fields->dz, xp->Dz2Dx2By[ix], 2, zDerivOrder);

  /* Dz^3 (Dx^2 By) */
  xp->Dz3Dx2By = (double**)czarray_2d(sizeof(double), fields->nx, fields->nz);
  for (ix=0; ix<fields->nx; ix++)
    takeMultipleDerivative(xp->Dx2By[ix], fields->nz, fields->dz, xp->Dz3Dx2By[ix], 3, zDerivOrder);

  /* Dz^2 (Dx^3 By) */
  xp->Dz2Dx3By = (double**)czarray_2d(sizeof(double), fields->nx, fields->nz);
  for (ix=0; ix<fields->nx; ix++)
    takeMultipleDerivative(xp->Dx3By[ix], fields->nz, fields->dz, xp->Dz2Dx3By[ix], 2, zDerivOrder);

  /* Dz^4 (Dx By) */
  xp->Dz4DxBy = (double**)czarray_2d(sizeof(double), fields->nx, fields->nz);
  for (ix=0; ix<fields->nx; ix++)
    takeMultipleDerivative(xp->DxBy[ix], fields->nz, fields->dz, xp->Dz4DxBy[ix], 4, zDerivOrder);

  /* Dz (Dx^4 By) */
  xp->DzDx4By = (double**)czarray_2d(sizeof(double), fields->nx, fields->nz);
  for (ix=0; ix<fields->nx; ix++)
    takeMultipleDerivative(xp->Dx4By[ix], fields->nz, fields->dz, xp->DzDx4By[ix], 1, zDerivOrder);

}

#include "matlib.h"

void SmoothData(MIDPLANE_FIELDS *fields, long smoothPoints, long fitTerms, double fraction, short versus_x, char *output) 
{
  MATRIX *A, *At, *AtA, *Y, *C, *S;
  long i, j, ix, iz;
  double factor;
  double *buffer;

  m_alloc(&A, smoothPoints, fitTerms);
  m_alloc(&At, fitTerms, smoothPoints);
  m_alloc(&AtA, fitTerms, fitTerms);
  m_alloc(&S, fitTerms, smoothPoints);
  m_alloc(&Y, smoothPoints, 1);
  m_alloc(&C, fitTerms, 1);

  for (i=0; i<smoothPoints; i++) {
    factor = 1;
    for (j=0; j<fitTerms; j++) {
      A->a[i][j] = factor;
      factor *= (i-smoothPoints/2);
    }
  }

  if (!m_trans(At, A) || !m_mult(AtA, At, A) || !m_invert(AtA, AtA) || !m_mult(S, AtA, At))
    SDDS_Bomb("matrix manipulation problem");
  m_free(&A);
  m_free(&At);
  m_free(&AtA);

  if (versus_x) {
    long ix0, dix;
    buffer = malloc(sizeof(*buffer)*fields->nx);
    for (iz=0; iz<fields->nz; iz++) {
      memset(buffer, 0, sizeof(*buffer)*fields->nx);
      for (ix=0; ix<fields->nx; ix++) {
        ix0 = ix-smoothPoints/2;
        dix = 0;
        if (ix0<0) {
          dix = ix0;
          ix0 = 0;
        }
        if ((ix0+smoothPoints-1)>=fields->nx) {
          ix0 = fields->nx - smoothPoints;
          dix = ix - (ix0 + smoothPoints/2);
        }
        for (i=0; i<smoothPoints; i++)
          Y->a[i][0] = fields->By[ix0+i][iz];
        m_mult(C, S, Y);
        if (dix==0) {
          /* usually true */
          buffer[ix] = C->a[0][0];
        } else {
          buffer[ix] = 0;
          factor = 1;
          for (i=0; i<fitTerms; i++) {
            buffer[ix] += C->a[i][0]*factor;
            factor *= dix;
          }
        }
      }

      for (ix=0; ix<fields->nx; ix++)
        fields->By[ix][iz] = fraction*buffer[ix] + (1-fraction)*fields->By[ix][iz];
    }
    free(buffer);
  } else {
    long iz0, diz;
    buffer = malloc(sizeof(*buffer)*fields->nz);
    for (ix=0; ix<fields->nx; ix++) {
      memset(buffer, 0, sizeof(*buffer)*fields->nz);
      for (iz=0; iz<fields->nz; iz++) {
        iz0 = iz-smoothPoints/2;
        diz = 0;
        if (iz0<0) {
          diz = iz0;
          iz0 = 0;
        }
        if ((iz0+smoothPoints-1)>=fields->nz) {
          iz0 = fields->nz - smoothPoints;
          diz = iz - (iz0 + smoothPoints/2);
        }
        for (i=0; i<smoothPoints; i++)
          Y->a[i][0] = fields->By[ix][iz0+i];
        m_mult(C, S, Y);
        if (diz==0) {
          /* usually true */
          buffer[iz] = C->a[0][0];
        } else {
          buffer[iz] = 0;
          factor = 1;
          for (i=0; i<fitTerms; i++) {
            buffer[iz] += C->a[i][0]*factor;
            factor *= diz;
          }
        }
      }

      for (iz=0; iz<fields->nz; iz++)
        fields->By[ix][iz] = fraction*buffer[iz] + (1-fraction)*fields->By[ix][iz];
    }
    free(buffer);
  }

  m_free(&S);
  m_free(&Y);
  m_free(&C);

  if (output) {
    SDDS_DATASET SDDSout;
    if (!SDDS_InitializeOutput(&SDDSout, SDDS_BINARY, 1, NULL, NULL, output) ||
        !SDDS_DefineSimpleParameter(&SDDSout, "x", "m", SDDS_DOUBLE) || 
        !SDDS_DefineSimpleColumn(&SDDSout, "z", "m", SDDS_DOUBLE) || 
        !SDDS_DefineSimpleColumn(&SDDSout, "By", "T", SDDS_DOUBLE) ||
        !SDDS_WriteLayout(&SDDSout))
      SDDS_Bomb("Problem writing post-smoothing data");

    for (ix=0; ix<fields->nx; ix++) {
      if (!SDDS_StartPage(&SDDSout, fields->nz) ||
          !SDDS_SetParameters(&SDDSout, SDDS_SET_BY_NAME|SDDS_PASS_BY_VALUE, "x", fields->x[ix], NULL) ||
          !SDDS_SetColumn(&SDDSout, SDDS_SET_BY_NAME, fields->z, fields->nz, "z") ||
          !SDDS_SetColumn(&SDDSout, SDDS_SET_BY_NAME, fields->By[ix], fields->nz, "By") ||
          !SDDS_WritePage(&SDDSout))
        SDDS_Bomb("Problem writing post-smoothing data");
    }
    if (!SDDS_Terminate(&SDDSout))
      SDDS_Bomb("Problem writing post-smoothing data");
  }
}
