
#include "mdb.h"
#include "SDDS.h"
#include "scan.h"
#include "fftpackC.h"
// #include "gsl/gsl_sf_bessel.h"

#define C_LIGHT_MKS 299792458.
#define E_CHARGE_PC 1.602176487e-7
#define E_MASS_MKS 9.1093837015e-31

typedef struct FRINGE_INT2 {
  double int1;
  double int2;
} FRINGE_INT2;
typedef struct FRINGE_INT3 {
  double int1;
  double int2;
  double int3;
} FRINGE_INT3;

typedef struct GGE_m_ell {
  double *grad;
} GGE_m_ell;

typedef struct COMPLEX {
  double re;
  double im;
} COMPLEX;

void readInGGE(char *ggeFile, double **z, double **ggeD, double **ggeQ, double **ggeS, int *Nz);
void readInAndMoveGGE(char *ggeFile, double **z, double **ggeD, double **ggeQ, double **ggeS, int *Nz, double xEntry);
void computeBOnCylinder(double xCenter, double radius, double **Bnorm, int *m,
                        GGE_m_ell **gges, int mMax, int l_2Max, int Nz, int Nangle);
void computeRequiredGGEs(double *dipole, double *DdipoleD2, double *quadrupole, double *sextupole,
                         double **Bnorm, double rho, double dz, int Nangle, int Nz);
void FFT(COMPLEX *field, int32_t isign, int32_t npts);
double Imp(double kR, long m);
double BesIn(double x, long order);
double BesselFuncIn(double x, double order);

void readInTrajectoryInput(char *input, double **x, double **xp, double **By, double *dz, int64_t *rows,
                           double *pCentral, int *Nsegments, double **zMagnetRef);

void LGBENDfringeCalc(double *z, double *ggeD, double *ggeQ, double *ggeS, double *xp,
                      double *x, double *By, int Nz, double invRigidity, double *zRef, int Nedges, char *output,
                      double totalBendAngle, char *edgeOutputFile, double zOffset);

void CCBENDfringeCalc(double *z, double *ggeD, double *ggeQ, double *ggeS,
                      int Nz, double invRigidity, double zRef, double bendAngle, double arcLength,
                      char *output, char *elementName);

// int refineNextMaximum(double *ggeD, int *zMaxInt, int edgeNum);

double integrateTrap(double *integrand, int startPt, int endPt, double dz);

void setStepFunction(double *heavside, double *maxD, int *zMaxInt, int *zEdgeInt, int edgeNum, double checkInt, double dz);

double intPower(double x, int n);

FRINGE_INT3 computeDipoleFringeInt(double *ggeD, double *heavside, double *z, double *zEdge,
                                   int *zMaxInt, int edgeNum, double invRigidity, double dz);

FRINGE_INT3 computeQuadrupoleFringeInt(double *ggeQ, double *heavside, double *z, double *zEdge,
                                       int *zMaxInt, int edgeNum, double invRigidity, double dz);

FRINGE_INT3 computeSextupoleFringeInt(double *ggeS, double *stepFuncS, double *z, double *zEdge,
                                      int *zMaxInt, int edgeNum, double invRigidity, double dz);

FRINGE_INT3 computeDipSextFringeInt(double *ggeD, double *stepFuncD, double *ggeS, double *stepFuncS,
                                    double *z, double *zEdge, int *zMaxInt, int edgeNum, double invRigidity, double dz);

#define SET_GGE 0
#define SET_CCBEND 1
#define SET_LGBEND 2
#define SET_PCENTRAL 3
#define SET_BENDANGLE 4
#define SET_ZREFFIELD 5
#define SET_TRAJECTORY 6
#define SET_ELEMENT_NAME 7
#define N_OPTIONS 8
char *option[N_OPTIONS] = {
  "gge",
  "ccbend",
  "lgbend",
  "pCentral",
  "bendAngle",
  "zRefField",
  "trajectory",
  "elementName",
};

#define USAGE "straightDipoleFringeCalc -gge=<filename> <outputfile>\n\
   -elementName=<string>\n\
   { -ccbend=pCentral=<value>,bendAngle=<radians>,[arcLength=<meters>,]xEntry=<meters>,zRefField=<meters>\n\
      | -lgbend=trajectory=<filename>,bendAngle=<radians>[,edgeOutput=<filename>]\n\
    }\n\
\n\
-gge         Give generalized gradient expansion file, such as produced by computeRBGGE\n\
             or computeCBGGE.\n\
<outputFile> Name of file to which to write fringe parameters in the format expected by\n\
             elegant's load_parameters command\n\
-elementName Name of the element to which the data will be loaded in elegant.\n\
-ccbend      CCBEND mode. This requires four additional parameters: pCentral, the\n\
             reference momentum (beta*gamma); bendAngle (with sign), in radians;\n\
             xEntry, the x coordinate of the hard edge magnet entry (and exit), along which\n\
             the fringe integrals will be calculated (defaults to 0); zRefField, the z location\n\
             in gge file from which the reference magnetic field values will be taken\n\
             (defaults to 0).\n\
-lgbend      LGBEND mode, which requires one additional parameter: \n\
             trajectory, name of an SDDS file from BGGEXP's PARTICLE_OUTPUT_FILE feature,\n\
             from tracking a single reference particle. The file should have some additional\n\
             parameters, which must be added by the user:\n\
             nMagnetSegments   Number of flat steps of non-zero field\n\
             zReference#       z location in file where reference\n\
                               magnetic field values will be taken,\n\
                               1<= # <= nMagnetSegments\n"

#define PCENTRAL_GIVEN 0x001UL
#define BENDANGLE_GIVEN 0x002UL
#define ZREFFIELD_GIVEN 0x004UL
#define TRAJECTORY_GIVEN 0x008UL
#define XENTRY_GIVEN 0x010UL
#define ARCLENGTH_GIVEN 0x020UL

#define ZREFMAX 1024

#define CCBEND_MODE 1
#define LGBEND_MODE 2

double zVertex = 0, zEntry = 0, zExit = 0;
double xVertex = 0, xEntry = 0, xExit = 0;

int main(int argc, char **argv) {
  SCANNED_ARG *scanned;
  long i_arg;
  unsigned long ccbendFlags, lgbendFlags;

  double *z = NULL, *ggeD = NULL, *ggeQ = NULL, *ggeS = NULL, *xp, *x, *By;
  double *zMagnetRef;

  double pCentral = DBL_MAX, invRigidity, bendAngle = DBL_MAX, arcLength = DBL_MAX;
  // 4.501990914651359e-05 - 1.922200520833333e-05; // should default to zero and be read in below

  int Nsegments, Nedges, ip, Nz;

  char *ggeFile = NULL, *output = NULL, *trajectoryFile = NULL, *elementName = NULL;
  char *edgeOutputFile = NULL;
  int mode = 0;

  SDDS_RegisterProgramName(argv[0]);

  zMagnetRef = malloc(sizeof(double) * ZREFMAX);
  for (ip = 0; ip < ZREFMAX; ip++)
    zMagnetRef[ip] = DBL_MAX;

  argc = scanargs(&scanned, argc, argv);
  if (argc < 2 || argc > (2 + N_OPTIONS)) {
    fprintf(stderr, "%s\n", USAGE);
    return (1);
  }
  for (i_arg = 1; i_arg < argc; i_arg++) {
    if (scanned[i_arg].arg_type == OPTION) {
      /* process options here */
      switch (match_string(scanned[i_arg].list[0], option, N_OPTIONS, 0)) {
      case SET_GGE:
        if (scanned[i_arg].n_items != 2) {
          fprintf(stderr, "invalid -gge syntax\n%s\n", USAGE);
          return (1);
        }
        ggeFile = scanned[i_arg].list[1];
        break;
      case SET_CCBEND:
        if (mode != 0)
          bomb("too many modes given", USAGE);
        mode = CCBEND_MODE;
        scanned[i_arg].n_items--;
        zMagnetRef[0] = 0;
        if (!scanItemList(&ccbendFlags, scanned[i_arg].list + 1, &scanned[i_arg].n_items, 0,
                          "pcentral", SDDS_DOUBLE, &pCentral, 1, PCENTRAL_GIVEN,
                          "bendangle", SDDS_DOUBLE, &bendAngle, 1, BENDANGLE_GIVEN,
                          "arclength", SDDS_DOUBLE, &arcLength, 1, ARCLENGTH_GIVEN,
                          "zreffield", SDDS_DOUBLE, &(zMagnetRef[0]), 1, ZREFFIELD_GIVEN,
                          "xentry", SDDS_DOUBLE, &(xEntry), 1, XENTRY_GIVEN,
                          NULL))
          bomb("invalid -ccbend syntax", USAGE);
        if (!(ccbendFlags & PCENTRAL_GIVEN))
          bomb("give pCentral value for -ccbend", USAGE);
        if (!(ccbendFlags & BENDANGLE_GIVEN))
          bomb("give bend angle value for -ccbend", USAGE);
        break;
      case SET_LGBEND:
        if (mode != 0)
          bomb("too many modes given", USAGE);
        mode = LGBEND_MODE;
        trajectoryFile = NULL;
        edgeOutputFile = NULL;
        scanned[i_arg].n_items--;
        lgbendFlags = 0;
        if (!scanItemList(&lgbendFlags, scanned[i_arg].list + 1, &scanned[i_arg].n_items, 0,
                          "trajectory", SDDS_STRING, &trajectoryFile, 1, TRAJECTORY_GIVEN,
                          "bendangle", SDDS_DOUBLE, &bendAngle, 1, BENDANGLE_GIVEN,
                          "edgeoutput", SDDS_STRING, &edgeOutputFile, 1, 0,
                          NULL) ||
            !(lgbendFlags & TRAJECTORY_GIVEN) || !trajectoryFile || !strlen(trajectoryFile) || !(lgbendFlags & BENDANGLE_GIVEN))
          bomb("invalid -lgbend syntax", USAGE);
        break;
      case SET_ELEMENT_NAME:
        if (scanned[i_arg].n_items != 2)
          bomb("invalid -elementName syntax", USAGE);
        elementName = scanned[i_arg].list[1];
        break;
      default:
        fprintf(stderr, "unknown option given\n%s\n", USAGE);
        return (1);
        break;
      }
    } else {
      if (!output)
        output = scanned[i_arg].list[0];
      else {
        fprintf(stderr, "too many filenames given\n%s\n", USAGE);
        return (1);
      }
    }
  }
  if (elementName == NULL)
    bomb("Must give -elementName option", NULL);
  if (ggeFile == NULL)
    bomb("Must include -gge option", NULL);
  if (output == NULL)
    bomb("No output file given", NULL);
  if (mode == 0)
    bomb("Must give -ccbend or -lgbend option", NULL);

  if (mode == CCBEND_MODE) {
    /* get GGE input */
    if (fabs(xEntry) > 1.0e-10) // read in and compute GGEs about new x=0
      readInAndMoveGGE(ggeFile, &z, &ggeD, &ggeQ, &ggeS, &Nz, xEntry);
    else
      readInGGE(ggeFile, &z, &ggeD, &ggeQ, &ggeS, &Nz);
    /****** Note that ggeD = -CnmS0[m=1]; ggeQ = -2*CnmS0[m=2]; ggeS = -6*CnmS0[m=3] + 0.25*CnmS2[m=1] ******/
    Nsegments = 1;
    Nedges = Nsegments + 1;
    zMagnetRef[0] -= z[0];
    for (ip = Nz - 1; ip >= 0; ip--)
      z[ip] -= z[0];
    /* add 0.01% of dz to eliminate rounding errors when converting to ints */
    // zMagnetRef[0] += 1.0e-4 * (z[1] - z[0]);
    invRigidity = 1.0e-12 * E_CHARGE_PC / (E_MASS_MKS * C_LIGHT_MKS * pCentral);
    CCBENDfringeCalc(z, ggeD, ggeQ, ggeS, Nz, invRigidity, zMagnetRef[0],
                     bendAngle, ccbendFlags & ARCLENGTH_GIVEN ? arcLength : -1, output, elementName);
  } else if (mode == LGBEND_MODE) {
    int64_t trajRows, i;
    double dzTraj, zOffset;

    readInGGE(ggeFile, &z, &ggeD, &ggeQ, &ggeS, &Nz);

    if (!trajectoryFile || !strlen(trajectoryFile))
      bomb("trajectory file must be given in LGBEND mode", USAGE);

    readInTrajectoryInput(trajectoryFile, &x, &xp, &By, &dzTraj, &trajRows, &pCentral, &Nsegments, &zMagnetRef);
    if (trajRows != (Nz + 1)) {
      fprintf(stderr, "Error: %ld rows in file %s but %ld rows in file %s\n",
              (long)trajRows, trajectoryFile, (long)Nz, ggeFile);
      exit(1);
    }
    for (i = 0; i < Nz - 1; i++)
      if (fabs(dzTraj - (z[i + 1] - z[i])) > 1e-6 * dzTraj) {
        fprintf(stderr, "dz values are mismatched between the GGE and trajectory files\n");
        exit(1);
      }

    Nedges = Nsegments + 1;
    printf("Nedges = %i\n", Nedges);
    /* Internally set z[0] = 0 to simplify calculations */
    zOffset = z[0];
    for (ip = 0; ip < Nedges - 1; ip++)
      zMagnetRef[ip] -= zOffset;
    for (ip = Nz - 1; ip >= 0; ip--)
      z[ip] -= zOffset;
    /* add 0.01% of dz to eliminate rounding errors when converting to ints */
    for (ip = 0; ip < Nedges - 1; ip++)
      zMagnetRef[ip] += 1.0e-4 * (z[1] - z[0]);
    invRigidity = 1.0e-12 * E_CHARGE_PC / (E_MASS_MKS * C_LIGHT_MKS * pCentral);
    zEntry -= zOffset;
    zVertex -= zOffset;
    zExit -= zOffset;
    LGBENDfringeCalc(z, ggeD, ggeQ, ggeS, xp, x, By, Nz, invRigidity, zMagnetRef, Nedges, output, bendAngle, edgeOutputFile,
                     zOffset);
  } else
    bomb("Something impossible happened.", NULL);

  // fclose(outFP);
  free(z);
  free(ggeD);
  free(ggeQ);
  free(ggeS);
  exit(0);
}

void LGBENDfringeCalc(double *z, double *ggeD, double *ggeQ, double *ggeS, double *xp,
                      double *x, double *By, int Nz, double invRigidity, double *zRef, int Nedges, char *output,
                      double totalBendAngle, char *edgeOutputFile, double zOffset) {
  SDDS_DATASET SDDSout;
  FILE *fpEdge = NULL;
  double *stepFuncD, *stepFuncQ, *stepFuncS;
  double *zEdge, *maxD, *maxQ, *maxS;
  double *bendRad, *edgeAngle, *edgeX;

  FRINGE_INT3 dipFringeInt, sextFringeInt;
  FRINGE_INT3 quadFringeInt, dipSextFringeInt;

  double temp1, fringeInt, dz, bendAngle;
  double integratedDipole;

  int *zEdgeInt, *zMaxInt;
  int edgeNum, ip, iRow;
  int NsegmentParams = 4;
  int NfringeParams = 10;

  if (edgeOutputFile) {
    fpEdge = fopen(edgeOutputFile, "w");
    fprintf(fpEdge, "SDDS1\n&column name=z type=double units=m &end\n");
    fprintf(fpEdge, "&column name=x type=double units=m &end\n");
    fprintf(fpEdge, "&column name=By type=double units=T &end\n");
    fprintf(fpEdge, "&data mode=ascii no_row_counts=1 &end\n");
  }

  dz = z[1] - z[0];
  /* allocate memory */
  stepFuncD = calloc(Nz, sizeof(double));
  stepFuncQ = calloc(Nz, sizeof(double));
  stepFuncS = calloc(Nz, sizeof(double));

  /* allocate memory */
  maxD = calloc(Nedges + 1, sizeof(double));
  maxQ = calloc(Nedges + 1, sizeof(double));
  maxS = calloc(Nedges + 1, sizeof(double));
  zEdge = calloc(Nedges, sizeof(double));
  bendRad = calloc(Nedges, sizeof(double));
  edgeAngle = calloc(Nedges, sizeof(double));
  edgeX = calloc(Nedges, sizeof(double));
  zEdgeInt = calloc(Nedges, sizeof(int));
  zMaxInt = calloc(Nedges + 1, sizeof(int));

  zMaxInt[0] = 0;
  for (edgeNum = 1; edgeNum < Nedges; edgeNum++) {
    zMaxInt[edgeNum] = (int)(zRef[edgeNum - 1] / dz);
    if ((zMaxInt[edgeNum] < 1) || (zMaxInt[edgeNum] > Nz - 2)) {
      printf("Magnetic reference point is located outside the field map!\n");
      return;
    }
  }
  zMaxInt[Nedges] = Nz - 1;

  /* define multipole content for each body segment including outside magnet */
  for (edgeNum = 0; edgeNum < Nedges + 1; edgeNum++) {
    maxD[edgeNum] = ggeD[zMaxInt[edgeNum]];
    maxQ[edgeNum] = ggeQ[zMaxInt[edgeNum]];
    maxS[edgeNum] = ggeS[zMaxInt[edgeNum]];
    /* set multipole content outside magnet to small value */
    if ((edgeNum == 0) || (edgeNum == Nedges)) {
      maxD[edgeNum] = 1.0e-10 * maxD[edgeNum];
      maxQ[edgeNum] = 1.0e-10 * maxQ[edgeNum];
      maxS[edgeNum] = 1.0e-10 * maxS[edgeNum];
    }
  }

  if (!SDDS_InitializeOutput(&SDDSout, SDDS_ASCII, 1, NULL, "Fringe integrals and body parameters for LGBEND", output)) {
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
    exit(1);
  }
  if (!SDDS_DefineSimpleColumn(&SDDSout, "ParameterName", NULL, SDDS_STRING) ||
      !SDDS_DefineSimpleColumn(&SDDSout, "ParameterValue", NULL, SDDS_DOUBLE) ||
      !SDDS_WriteLayout(&SDDSout)) {
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
    exit(1);
  }

  edgeNum = 0;

  /* compute the integrated bending field */
  integratedDipole = integrateTrap(ggeD, 0, Nz - 1, dz);
  /* find hard edge according to integrated Bfield */
  fringeInt = integrateTrap(ggeD, zMaxInt[edgeNum], zMaxInt[edgeNum + 1], dz);
  zEdge[edgeNum] = (maxD[edgeNum + 1] * (z[zMaxInt[edgeNum + 1]] - z[zMaxInt[edgeNum]]) - fringeInt) / (maxD[edgeNum + 1] - maxD[edgeNum]);
  zEdge[edgeNum] += z[zMaxInt[edgeNum]];
  temp1 = (zEdge[edgeNum] - z[(int)(zEdge[edgeNum] / dz)]) / dz;
  /* angle at first edge set by input x' from trajectory */
  edgeAngle[edgeNum] = atan(xp[0]);
  /* horizontal coordinate set by linear interpolation to edge */
  edgeX[edgeNum] = (1.0 - temp1) * x[(int)(zEdge[edgeNum] / dz)] + temp1 * x[(int)(zEdge[edgeNum] / dz) + 1];
  printf("zEdge[%d] = %le, edgeX[%d] = %le\n", edgeNum, zEdge[edgeNum], edgeNum, edgeX[edgeNum]);
  if (fpEdge)
    fprintf(fpEdge, "%le %le %le\n", zEdge[edgeNum] + zOffset, edgeX[edgeNum], By[(int)(zEdge[edgeNum] / dz)]);

  /* initial coordinates (to be set by command line input) */
  /*
    edgeAngle[edgeNum] = atan(xp[0]);
    edgeX[edgeNum] = -3.003057073789646e-03;
    printf("edgeX[%d] = %le\n", edgeNum, edgeX[edgeNum]);
  */

  bendRad[edgeNum] = invRigidity * maxD[edgeNum + 1];
  bendRad[edgeNum] = 1.0 / bendRad[edgeNum];

  /* find step function field that gives same integrated field */
  zEdgeInt[edgeNum] = (int)(zEdge[edgeNum] / dz);
  setStepFunction(stepFuncD, maxD, zMaxInt, zEdgeInt, edgeNum, fringeInt, dz);

  for (ip = zMaxInt[edgeNum]; ip < zEdgeInt[edgeNum]; ip++)
    stepFuncQ[ip] = maxQ[edgeNum];
  stepFuncQ[zEdgeInt[edgeNum]] = ((maxQ[edgeNum] - maxQ[edgeNum + 1]) * stepFuncD[zEdgeInt[edgeNum]] - maxQ[edgeNum] * maxD[edgeNum + 1] + maxQ[edgeNum + 1] * maxD[edgeNum]) / (maxD[edgeNum] - maxD[edgeNum + 1]);
  for (ip = zEdgeInt[edgeNum] + 1; ip <= zMaxInt[edgeNum + 1]; ip++)
    stepFuncQ[ip] = maxQ[edgeNum + 1];

  for (ip = zMaxInt[edgeNum]; ip < zEdgeInt[edgeNum]; ip++)
    stepFuncS[ip] = maxS[edgeNum];
  stepFuncS[zEdgeInt[edgeNum]] = ((maxS[edgeNum] - maxS[edgeNum + 1]) * stepFuncD[zEdgeInt[edgeNum]] - maxS[edgeNum] * maxD[edgeNum + 1] + maxS[edgeNum + 1] * maxD[edgeNum]) / (maxD[edgeNum] - maxD[edgeNum + 1]);
  for (ip = zEdgeInt[edgeNum] + 1; ip <= zMaxInt[edgeNum + 1]; ip++)
    stepFuncS[ip] = maxS[edgeNum + 1];

  dipFringeInt = computeDipoleFringeInt(ggeD, stepFuncD, z, zEdge, zMaxInt, edgeNum, invRigidity, dz);
  quadFringeInt = computeQuadrupoleFringeInt(ggeQ, stepFuncQ, z, zEdge, zMaxInt, edgeNum, invRigidity, dz);
  sextFringeInt = computeSextupoleFringeInt(ggeS, stepFuncS, z, zEdge, zMaxInt, edgeNum, invRigidity, dz);
  dipSextFringeInt = computeDipSextFringeInt(ggeD, stepFuncD, ggeS, stepFuncS, z, zEdge, zMaxInt, edgeNum, invRigidity, dz);

  // Debugging
  // printf("%23.15e %23.15e \n", zEdge[edgeNum], stepFuncD[zMaxInt[edgeNum]]);
  // printf("%23.15e %23.15e \n", zEdge[edgeNum], stepFuncD[zMaxInt[edgeNum + 1]]);
  // printf("%23.15e \n", edgeAngle[edgeNum]);

  for (edgeNum = 1; edgeNum < Nedges; edgeNum++) {
    /* find hard edge according to integrated Bfield */
    fringeInt = integrateTrap(ggeD, zMaxInt[edgeNum], zMaxInt[edgeNum + 1], dz);
    zEdge[edgeNum] = (maxD[edgeNum + 1] * (z[zMaxInt[edgeNum + 1]] - z[zMaxInt[edgeNum]]) - fringeInt) / (maxD[edgeNum + 1] - maxD[edgeNum]);
    zEdge[edgeNum] += z[zMaxInt[edgeNum]];
    temp1 = (zEdge[edgeNum] - z[(int)(zEdge[edgeNum] / dz)]) / dz;
    /* angles at last edge set by output x' from trajectory */
    if (edgeNum == Nedges - 1)
      edgeAngle[edgeNum] = atan(xp[Nz - 1]);
    else
      edgeAngle[edgeNum] = atan((1.0 - temp1) * xp[(int)(zEdge[edgeNum] / dz)] + temp1 * xp[(int)(zEdge[edgeNum] / dz) + 1]);
    /* interior angles set by linear interpolation to edge */
    /* horizontal coordinate set by linear interpolation to edge */
    edgeX[edgeNum] = (1.0 - temp1) * x[(int)(zEdge[edgeNum] / dz)] + temp1 * x[(int)(zEdge[edgeNum] / dz) + 1];

    bendRad[edgeNum] = invRigidity * maxD[edgeNum + 1];
    bendRad[edgeNum] = 1.0 / bendRad[edgeNum];

    /* Compute the fraction of the integrated bending field in this segment */
    bendAngle = maxD[edgeNum] * (zEdge[edgeNum] - zEdge[edgeNum - 1]) / integratedDipole;
    bendAngle = bendAngle * totalBendAngle;
    /* Compute next edge (x,x') (will be done by elegant) */
    edgeAngle[edgeNum] = edgeAngle[edgeNum - 1] - bendAngle;
    edgeX[edgeNum] = edgeX[edgeNum - 1] + (zEdge[edgeNum] - zEdge[edgeNum - 1]) * ((cos(edgeAngle[edgeNum]) - cos(edgeAngle[edgeNum - 1])) / (sin(edgeAngle[edgeNum - 1]) - sin(edgeAngle[edgeNum])));

    /* find step function field that gives same integrated field */
    zEdgeInt[edgeNum] = (int)(zEdge[edgeNum] / dz);
    setStepFunction(stepFuncD, maxD, zMaxInt, zEdgeInt, edgeNum, fringeInt, dz);

    if (fpEdge)
      fprintf(fpEdge, "%le %le %le\n", zEdge[edgeNum] + zOffset, edgeX[edgeNum], By[zEdgeInt[edgeNum]]);

    // Debugging
    // printf("%23.15e %23.15e \n", zEdge[edgeNum], stepFuncD[zMaxInt[edgeNum]]);
    // printf("%23.15e %23.15e \n", zEdge[edgeNum], stepFuncD[zMaxInt[edgeNum + 1]]);
    // printf("%23.15e \n", edgeAngle[edgeNum]);

    for (ip = zMaxInt[edgeNum]; ip < zEdgeInt[edgeNum]; ip++)
      stepFuncQ[ip] = maxQ[edgeNum];
    stepFuncQ[zEdgeInt[edgeNum]] = ((maxQ[edgeNum] - maxQ[edgeNum + 1]) * stepFuncD[zEdgeInt[edgeNum]] - maxQ[edgeNum] * maxD[edgeNum + 1] + maxQ[edgeNum + 1] * maxD[edgeNum]) / (maxD[edgeNum] - maxD[edgeNum + 1]);
    for (ip = zEdgeInt[edgeNum] + 1; ip <= zMaxInt[edgeNum + 1]; ip++)
      stepFuncQ[ip] = maxQ[edgeNum + 1];

    for (ip = zMaxInt[edgeNum]; ip < zEdgeInt[edgeNum]; ip++)
      stepFuncS[ip] = maxS[edgeNum];
    stepFuncS[zEdgeInt[edgeNum]] = ((maxS[edgeNum] - maxS[edgeNum + 1]) * stepFuncD[zEdgeInt[edgeNum]] - maxS[edgeNum] * maxD[edgeNum + 1] + maxS[edgeNum + 1] * maxD[edgeNum]) / (maxD[edgeNum] - maxD[edgeNum + 1]);
    for (ip = zEdgeInt[edgeNum] + 1; ip <= zMaxInt[edgeNum + 1]; ip++)
      stepFuncS[ip] = maxS[edgeNum + 1];

    iRow = 0;
    /* Output longitudinal length, bend angle, and multipole content */
    if (!SDDS_StartPage(&SDDSout, NsegmentParams + (edgeNum == 1 ? 2 : 1) * NfringeParams + 2 + (edgeNum == 1 ? 6 : 0)) ||
        !SDDS_SetRowValues(&SDDSout, SDDS_SET_BY_NAME | SDDS_PASS_BY_VALUE, iRow++,
                           "ParameterName", "LONGIT_L",
                           "ParameterValue", zEdge[edgeNum] - zEdge[edgeNum - 1], NULL) ||
        !SDDS_SetRowValues(&SDDSout, SDDS_SET_BY_NAME | SDDS_PASS_BY_VALUE, iRow++,
                           "ParameterName", "ANGLE",
                           "ParameterValue", bendAngle, NULL) ||
        //"ParameterValue", edgeAngle[edgeNum-1]-edgeAngle[edgeNum], NULL) ||
        !SDDS_SetRowValues(&SDDSout, SDDS_SET_BY_NAME | SDDS_PASS_BY_VALUE, iRow++,
                           "ParameterName", "K1",
                           "ParameterValue", invRigidity * maxQ[edgeNum], NULL) ||
        !SDDS_SetRowValues(&SDDSout, SDDS_SET_BY_NAME | SDDS_PASS_BY_VALUE, iRow++,
                           "ParameterName", "K2",
                           "ParameterValue", invRigidity * (maxS[edgeNum] - 0.25 * maxD[edgeNum]), NULL)) {
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
      exit(1);
    }

    /* Output the first edge x coordinate and fringe field terms if appropriate */
    if (edgeNum == 1) {
      if (!SDDS_SetRowValues(&SDDSout, SDDS_SET_BY_NAME | SDDS_PASS_BY_VALUE, iRow++,
                             "ParameterName", "ENTRY_X",
                             "ParameterValue", edgeX[edgeNum - 1], NULL) ||
          !SDDS_SetRowValues(&SDDSout, SDDS_SET_BY_NAME | SDDS_PASS_BY_VALUE, iRow++,
                             "ParameterName", "ENTRY_ANGLE",
                             "ParameterValue", edgeAngle[edgeNum - 1], NULL) ||
          !SDDS_SetRowValues(&SDDSout, SDDS_SET_BY_NAME | SDDS_PASS_BY_VALUE, iRow++,
                             "ParameterName", "FRINGE1K0",
                             "ParameterValue", dipFringeInt.int2, NULL) ||
          !SDDS_SetRowValues(&SDDSout, SDDS_SET_BY_NAME | SDDS_PASS_BY_VALUE, iRow++,
                             "ParameterName", "FRINGE1K2",
                             "ParameterValue", dipFringeInt.int3, NULL) ||
          !SDDS_SetRowValues(&SDDSout, SDDS_SET_BY_NAME | SDDS_PASS_BY_VALUE, iRow++,
                             "ParameterName", "FRINGE1K4",
                             "ParameterValue", sextFringeInt.int3, NULL) ||
          !SDDS_SetRowValues(&SDDSout, SDDS_SET_BY_NAME | SDDS_PASS_BY_VALUE, iRow++,
                             "ParameterName", "FRINGE1K5",
                             "ParameterValue", sextFringeInt.int2, NULL) ||
          !SDDS_SetRowValues(&SDDSout, SDDS_SET_BY_NAME | SDDS_PASS_BY_VALUE, iRow++,
                             "ParameterName", "FRINGE1K6",
                             "ParameterValue", sextFringeInt.int1, NULL) ||
          !SDDS_SetRowValues(&SDDSout, SDDS_SET_BY_NAME | SDDS_PASS_BY_VALUE, iRow++,
                             "ParameterName", "FRINGE1K7",
                             "ParameterValue", dipSextFringeInt.int1 + dipSextFringeInt.int2, NULL) ||
          !SDDS_SetRowValues(&SDDSout, SDDS_SET_BY_NAME | SDDS_PASS_BY_VALUE, iRow++,
                             "ParameterName", "FRINGE1I0",
                             "ParameterValue", quadFringeInt.int1, NULL) ||
          !SDDS_SetRowValues(&SDDSout, SDDS_SET_BY_NAME | SDDS_PASS_BY_VALUE, iRow++,
                             "ParameterName", "FRINGE1I1",
                             "ParameterValue", quadFringeInt.int2, NULL) ||
          !SDDS_SetRowValues(&SDDSout, SDDS_SET_BY_NAME | SDDS_PASS_BY_VALUE, iRow++,
                             "ParameterName", "PREDRIFT",
                             "ParameterValue", (zEdge[edgeNum - 1] - zEntry) / cos(edgeAngle[edgeNum - 1]), NULL) ||
          !SDDS_SetRowValues(&SDDSout, SDDS_SET_BY_NAME | SDDS_PASS_BY_VALUE, iRow++,
                             "ParameterName", "XEntry", "ParameterValue", xEntry, NULL) ||
          !SDDS_SetRowValues(&SDDSout, SDDS_SET_BY_NAME | SDDS_PASS_BY_VALUE, iRow++,
                             "ParameterName", "ZEntry", "ParameterValue", zEntry + zOffset, NULL) ||
          !SDDS_SetRowValues(&SDDSout, SDDS_SET_BY_NAME | SDDS_PASS_BY_VALUE, iRow++,
                             "ParameterName", "XVertex", "ParameterValue", xVertex, NULL) ||
          !SDDS_SetRowValues(&SDDSout, SDDS_SET_BY_NAME | SDDS_PASS_BY_VALUE, iRow++,
                             "ParameterName", "ZVertex", "ParameterValue", zVertex + zOffset, NULL) ||
          !SDDS_SetRowValues(&SDDSout, SDDS_SET_BY_NAME | SDDS_PASS_BY_VALUE, iRow++,
                             "ParameterName", "XExit", "ParameterValue", xExit, NULL) ||
          !SDDS_SetRowValues(&SDDSout, SDDS_SET_BY_NAME | SDDS_PASS_BY_VALUE, iRow++,
                             "ParameterName", "ZExit", "ParameterValue", zExit + zOffset, NULL)) {
        SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
        exit(1);
      }
    }

    dipFringeInt = computeDipoleFringeInt(ggeD, stepFuncD, z, zEdge, zMaxInt, edgeNum, invRigidity, dz);
    quadFringeInt = computeQuadrupoleFringeInt(ggeQ, stepFuncQ, z, zEdge, zMaxInt, edgeNum, invRigidity, dz);
    sextFringeInt = computeSextupoleFringeInt(ggeS, stepFuncS, z, zEdge, zMaxInt, edgeNum, invRigidity, dz);
    dipSextFringeInt = computeDipSextFringeInt(ggeD, stepFuncD, ggeS, stepFuncS, z, zEdge, zMaxInt, edgeNum, invRigidity, dz);

    // The terms EXIT_X and EXIT_ANGLE will be calculated by elegant, and are temporarily here for checking
    if (!SDDS_SetRowValues(&SDDSout, SDDS_SET_BY_NAME | SDDS_PASS_BY_VALUE, iRow++,
                           "ParameterName", "EXIT_X",
                           "ParameterValue", edgeX[edgeNum], NULL) ||
        !SDDS_SetRowValues(&SDDSout, SDDS_SET_BY_NAME | SDDS_PASS_BY_VALUE, iRow++,
                           "ParameterName", "EXIT_ANGLE",
                           "ParameterValue", edgeAngle[edgeNum], NULL) ||
        !SDDS_SetRowValues(&SDDSout, SDDS_SET_BY_NAME | SDDS_PASS_BY_VALUE, iRow++,
                           "ParameterName", "FRINGE2K0",
                           "ParameterValue", dipFringeInt.int2, NULL) ||
        !SDDS_SetRowValues(&SDDSout, SDDS_SET_BY_NAME | SDDS_PASS_BY_VALUE, iRow++,
                           "ParameterName", "FRINGE2K2",
                           "ParameterValue", dipFringeInt.int3, NULL) ||
        !SDDS_SetRowValues(&SDDSout, SDDS_SET_BY_NAME | SDDS_PASS_BY_VALUE, iRow++,
                           "ParameterName", "FRINGE2K4",
                           "ParameterValue", sextFringeInt.int3, NULL) ||
        !SDDS_SetRowValues(&SDDSout, SDDS_SET_BY_NAME | SDDS_PASS_BY_VALUE, iRow++,
                           "ParameterName", "FRINGE2K5",
                           "ParameterValue", sextFringeInt.int2, NULL) ||
        !SDDS_SetRowValues(&SDDSout, SDDS_SET_BY_NAME | SDDS_PASS_BY_VALUE, iRow++,
                           "ParameterName", "FRINGE2K6",
                           "ParameterValue", sextFringeInt.int1, NULL) ||
        !SDDS_SetRowValues(&SDDSout, SDDS_SET_BY_NAME | SDDS_PASS_BY_VALUE, iRow++,
                           "ParameterName", "FRINGE2K7",
                           "ParameterValue", dipSextFringeInt.int1 + dipSextFringeInt.int2, NULL) ||
        !SDDS_SetRowValues(&SDDSout, SDDS_SET_BY_NAME | SDDS_PASS_BY_VALUE, iRow++,
                           "ParameterName", "FRINGE2I0",
                           "ParameterValue", quadFringeInt.int1, NULL) ||
        !SDDS_SetRowValues(&SDDSout, SDDS_SET_BY_NAME | SDDS_PASS_BY_VALUE, iRow++,
                           "ParameterName", "FRINGE2I1",
                           "ParameterValue", quadFringeInt.int2, NULL) ||
        (edgeNum == (Nedges - 1) &&
         !SDDS_SetRowValues(&SDDSout, SDDS_SET_BY_NAME | SDDS_PASS_BY_VALUE, iRow++,
                            "ParameterName", "POSTDRIFT",
                            "ParameterValue", (zExit - zEdge[edgeNum]) / cos(edgeAngle[edgeNum]), NULL))) {
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
      exit(1);
    }
    if (!SDDS_WritePage(&SDDSout)) {
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
      exit(1);
    }
  }
  if (!SDDS_Terminate(&SDDSout)) {
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
    exit(1);
  }
  free(stepFuncD);
  free(stepFuncQ);
  free(stepFuncS);
  free(maxD);
  free(maxQ);
  free(maxS);
  free(zEdge);
  free(bendRad);
  free(edgeAngle);
  free(edgeX);
  free(zEdgeInt);
  free(zMaxInt);
}

void CCBENDfringeCalc(double *z, double *ggeD, double *ggeQ, double *ggeS, int Nz, double invRigidity,
                      double zRef, double bendAngle, double arcLength, char *output, char *elementName) {
  SDDS_DATASET SDDSout;
  double *stepFuncD, *stepFuncQ, *stepFuncS;
  double *zEdge, *maxD, *maxQ, *maxS;

  FRINGE_INT3 dipFringeInt, sextFringeInt;
  FRINGE_INT3 quadFringeInt, dipSextFringeInt;

  double fringeInt, dz = z[1] - z[0];

  int *zEdgeInt, *zMaxInt;
  int ip, edgeNum, Nedges = 2, iRow;

  /* allocate memory */
  stepFuncD = calloc(Nz, sizeof(double));
  stepFuncQ = calloc(Nz, sizeof(double));
  stepFuncS = calloc(Nz, sizeof(double));

  maxD = calloc(Nedges + 1, sizeof(double));
  maxQ = calloc(Nedges + 1, sizeof(double));
  maxS = calloc(Nedges + 1, sizeof(double));
  zEdge = calloc(Nedges, sizeof(double));

  zEdgeInt = calloc(Nedges, sizeof(int));
  zMaxInt = calloc(Nedges + 1, sizeof(int));

  /* set position of the reference points */
  zMaxInt[0] = 0;
  zMaxInt[1] = (int)(zRef / dz);
  if ((zMaxInt[1] < 1) || (zMaxInt[1] > Nz - 2)) {
    printf("Magnetic reference point is located outside the field map!\n");
    return;
  }
  zMaxInt[2] = Nz - 1;

  /* set multipole content outside magnet to small value */
  maxD[0] = 1.0e-10 * ggeD[zMaxInt[0]];
  maxQ[0] = 1.0e-10 * ggeQ[zMaxInt[0]];
  maxS[0] = 1.0e-10 * ggeS[zMaxInt[0]];
  maxD[2] = 1.0e-10 * ggeD[zMaxInt[2]];
  maxQ[2] = 1.0e-10 * ggeQ[zMaxInt[2]];
  maxS[2] = 1.0e-10 * ggeS[zMaxInt[2]];

  /* set multipole content in body */
  maxD[1] = ggeD[zMaxInt[1]];
  maxQ[1] = ggeQ[zMaxInt[1]];
  maxS[1] = ggeS[zMaxInt[1]];

  if (!SDDS_InitializeOutput(&SDDSout, SDDS_ASCII, 1, NULL, "Fringe integrals and body parameters for CCBEND", output)) {
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
    exit(1);
  }
  if (!SDDS_DefineSimpleColumn(&SDDSout, "ElementName", NULL, SDDS_STRING) ||
      !SDDS_DefineSimpleColumn(&SDDSout, "ElementType", NULL, SDDS_STRING) ||
      !SDDS_DefineSimpleColumn(&SDDSout, "ElementParameter", NULL, SDDS_STRING) ||
      !SDDS_DefineSimpleColumn(&SDDSout, "ParameterValue", NULL, SDDS_DOUBLE) ||
      !SDDS_WriteLayout(&SDDSout) ||
      !SDDS_StartPage(&SDDSout, 12 * Nedges)) {
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
    exit(1);
  }

  // Should find a better way to do this...
  fringeInt = integrateTrap(ggeD, zMaxInt[0], zMaxInt[1], dz);
  zEdge[0] = (maxD[1] * (z[zMaxInt[1]] - z[zMaxInt[0]]) - fringeInt) / (maxD[1] - maxD[0]);
  zEdge[0] += z[zMaxInt[0]];
  fringeInt = integrateTrap(ggeD, zMaxInt[1], zMaxInt[2], dz);
  zEdge[1] = (maxD[2] * (z[zMaxInt[2]] - z[zMaxInt[1]]) - fringeInt) / (maxD[2] - maxD[1]);
  zEdge[1] += z[zMaxInt[1]];
  if (arcLength < 0)
    /* compute the approximate arc length if the user didn't give it */
    arcLength = 0.5 * bendAngle * (zEdge[1] - zEdge[0]) / sin(0.5 * bendAngle);

  for (edgeNum = iRow = 0; edgeNum < Nedges; edgeNum++) {
    char nameBuffer[256];
    /* find hard edge according to integrated Bfield */
    fringeInt = integrateTrap(ggeD, zMaxInt[edgeNum], zMaxInt[edgeNum + 1], dz);
    zEdge[edgeNum] = (maxD[edgeNum + 1] * (z[zMaxInt[edgeNum + 1]] - z[zMaxInt[edgeNum]]) - fringeInt) / (maxD[edgeNum + 1] - maxD[edgeNum]);
    zEdge[edgeNum] += z[zMaxInt[edgeNum]];

    /* find step function field that gives same integrated field */
    zEdgeInt[edgeNum] = (int)(zEdge[edgeNum] / dz);
    setStepFunction(stepFuncD, maxD, zMaxInt, zEdgeInt, edgeNum, fringeInt, dz);

    for (ip = zMaxInt[edgeNum]; ip < zEdgeInt[edgeNum]; ip++)
      stepFuncQ[ip] = maxQ[edgeNum];
    stepFuncQ[zEdgeInt[edgeNum]] = ((maxQ[edgeNum] - maxQ[edgeNum + 1]) * stepFuncD[zEdgeInt[edgeNum]] - maxQ[edgeNum] * maxD[edgeNum + 1] + maxQ[edgeNum + 1] * maxD[edgeNum]) / (maxD[edgeNum] - maxD[edgeNum + 1]);
    for (ip = zEdgeInt[edgeNum] + 1; ip <= zMaxInt[edgeNum + 1]; ip++)
      stepFuncQ[ip] = maxQ[edgeNum + 1];

    for (ip = zMaxInt[edgeNum]; ip < zEdgeInt[edgeNum]; ip++)
      stepFuncS[ip] = maxS[edgeNum];
    stepFuncS[zEdgeInt[edgeNum]] = ((maxS[edgeNum] - maxS[edgeNum + 1]) * stepFuncD[zEdgeInt[edgeNum]] - maxS[edgeNum] * maxD[edgeNum + 1] + maxS[edgeNum + 1] * maxD[edgeNum]) / (maxD[edgeNum] - maxD[edgeNum + 1]);
    for (ip = zEdgeInt[edgeNum] + 1; ip <= zMaxInt[edgeNum + 1]; ip++)
      stepFuncS[ip] = maxS[edgeNum + 1];

    dipFringeInt = computeDipoleFringeInt(ggeD, stepFuncD, z, zEdge, zMaxInt, edgeNum, invRigidity, dz);
    quadFringeInt = computeQuadrupoleFringeInt(ggeQ, stepFuncQ, z, zEdge, zMaxInt, edgeNum, invRigidity, dz);
    sextFringeInt = computeSextupoleFringeInt(ggeS, stepFuncS, z, zEdge, zMaxInt, edgeNum, invRigidity, dz);

    dipSextFringeInt = computeDipSextFringeInt(ggeD, stepFuncD, ggeS, stepFuncS, z, zEdge, zMaxInt, edgeNum, invRigidity, dz);

    // The arc length is calculated above...
    // if (edgeNum == 1)
    // arcLength = 0.5*bendAngle*(zEdge[1]-zEdge[0])/sin(0.5*bendAngle);
    // else
    // arcLength = bendAngle / (invRigidity * maxD[1]);

    if (edgeNum == 0) {
      if (!SDDS_SetRowValues(&SDDSout, SDDS_SET_BY_NAME | SDDS_PASS_BY_VALUE, iRow++,
                             "ElementName", elementName, "ElementType", "CCBEND",
                             "ElementParameter", "L",
                             "ParameterValue", arcLength, NULL) ||
          !SDDS_SetRowValues(&SDDSout, SDDS_SET_BY_NAME | SDDS_PASS_BY_VALUE, iRow++,
                             "ElementName", elementName, "ElementType", "CCBEND",
                             "ElementParameter", "ANGLE",
                             "ParameterValue", bendAngle, NULL) ||
          !SDDS_SetRowValues(&SDDSout, SDDS_SET_BY_NAME | SDDS_PASS_BY_VALUE, iRow++,
                             "ElementName", elementName, "ElementType", "CCBEND",
                             "ElementParameter", "K1",
                             "ParameterValue", invRigidity * maxQ[edgeNum + 1], NULL) ||
          !SDDS_SetRowValues(&SDDSout, SDDS_SET_BY_NAME | SDDS_PASS_BY_VALUE, iRow++,
                             "ElementName", elementName, "ElementType", "CCBEND",
                             "ElementParameter", "K2",
                             "ParameterValue", invRigidity * (maxS[edgeNum + 1] - 0.25 * maxD[edgeNum + 1]), NULL) ||
          !SDDS_SetRowValues(&SDDSout, SDDS_SET_BY_NAME | SDDS_PASS_BY_VALUE, iRow++,
                             "ElementName", elementName, "ElementType", "CCBEND",
                             "ElementParameter", "FRINGEMODEL",
                             "ParameterValue", (double)1.0, NULL) ||
          !SDDS_SetRowValues(&SDDSout, SDDS_SET_BY_NAME | SDDS_PASS_BY_VALUE, iRow++,
                             "ElementName", elementName, "ElementType", "CCBEND",
                             "ElementParameter", "COMPENSATE_KN",
                             "ParameterValue", (double)1.0, NULL))

        {
          SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
          exit(1);
        }
    }

    if (!SDDS_SetRowValues(&SDDSout, SDDS_SET_BY_NAME | SDDS_PASS_BY_VALUE, iRow++,
                           "ElementName", elementName, "ElementType", "CCBEND",
                           "ElementParameter", snprintf(nameBuffer, 256, "FRINGE%dK0", edgeNum + 1) ? nameBuffer : "?",
                           "ParameterValue", dipFringeInt.int2, NULL) ||
        !SDDS_SetRowValues(&SDDSout, SDDS_SET_BY_NAME | SDDS_PASS_BY_VALUE, iRow++,
                           "ElementName", elementName, "ElementType", "CCBEND",
                           "ElementParameter", snprintf(nameBuffer, 256, "FRINGE%dK2", edgeNum + 1) ? nameBuffer : "?",
                           "ParameterValue", dipFringeInt.int3, NULL) ||
        !SDDS_SetRowValues(&SDDSout, SDDS_SET_BY_NAME | SDDS_PASS_BY_VALUE, iRow++,
                           "ElementName", elementName, "ElementType", "CCBEND",
                           "ElementParameter", snprintf(nameBuffer, 256, "FRINGE%dK4", edgeNum + 1) ? nameBuffer : "?",
                           "ParameterValue", sextFringeInt.int3, NULL) ||
        !SDDS_SetRowValues(&SDDSout, SDDS_SET_BY_NAME | SDDS_PASS_BY_VALUE, iRow++,
                           "ElementName", elementName, "ElementType", "CCBEND",
                           "ElementParameter", snprintf(nameBuffer, 256, "FRINGE%dK5", edgeNum + 1) ? nameBuffer : "?",
                           "ParameterValue", sextFringeInt.int2, NULL) ||
        !SDDS_SetRowValues(&SDDSout, SDDS_SET_BY_NAME | SDDS_PASS_BY_VALUE, iRow++,
                           "ElementName", elementName, "ElementType", "CCBEND",
                           "ElementParameter", snprintf(nameBuffer, 256, "FRINGE%dK6", edgeNum + 1) ? nameBuffer : "?",
                           "ParameterValue", sextFringeInt.int1, NULL) ||
        !SDDS_SetRowValues(&SDDSout, SDDS_SET_BY_NAME | SDDS_PASS_BY_VALUE, iRow++,
                           "ElementName", elementName, "ElementType", "CCBEND",
                           "ElementParameter", snprintf(nameBuffer, 256, "FRINGE%dK7", edgeNum + 1) ? nameBuffer : "?",
                           "ParameterValue", dipSextFringeInt.int1 + dipSextFringeInt.int2, NULL) ||
        !SDDS_SetRowValues(&SDDSout, SDDS_SET_BY_NAME | SDDS_PASS_BY_VALUE, iRow++,
                           "ElementName", elementName, "ElementType", "CCBEND",
                           "ElementParameter", snprintf(nameBuffer, 256, "FRINGE%dI0", edgeNum + 1) ? nameBuffer : "?",
                           "ParameterValue", quadFringeInt.int1, NULL) ||
        !SDDS_SetRowValues(&SDDSout, SDDS_SET_BY_NAME | SDDS_PASS_BY_VALUE, iRow++,
                           "ElementName", elementName, "ElementType", "CCBEND",
                           "ElementParameter", snprintf(nameBuffer, 256, "FRINGE%dI1", edgeNum + 1) ? nameBuffer : "?",
                           "ParameterValue", quadFringeInt.int2, NULL)) {
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
      exit(1);
    }
  }
  if (!SDDS_WritePage(&SDDSout) || !SDDS_Terminate(&SDDSout)) {
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
    exit(1);
  }
  printf("  zEdge1 = %23.15e, zEdge2 = %23.15e\n", zEdge[0], zEdge[1]);
  printf("Difference zEdge2 - zEdge2 = %23.15e\n", zEdge[1] - zEdge[0]);
  printf("                Arc Length = %23.15e\n", 0.5 * bendAngle * (zEdge[1] - zEdge[0]) / sin(0.5 * bendAngle));
  free(stepFuncD);
  free(stepFuncQ);
  free(stepFuncS);
  free(maxD);
  free(maxQ);
  free(maxS);
  free(zEdge);
  free(zEdgeInt);
  free(zMaxInt);
}

void readInGGE(char *ggeFile, double **z, double **ggeD, double **ggeQ, double **ggeS, int *Nz) {
  SDDS_DATASET SDDSin;
  int64_t i, rows;
  double *zValues, *CnmS0Values, *CnmS2Values;

  if (!SDDS_InitializeInput(&SDDSin, ggeFile)) {
    fprintf(stderr, "Unable to read file %s\n", ggeFile);
    exit(1);
  }
  if (SDDS_CheckColumn(&SDDSin, "z", "m", SDDS_ANY_FLOATING_TYPE, stderr) != SDDS_CHECK_OK) {
    fprintf(stderr, "Unable to find floating-point column \"z\" with units \"m\" in file %s\n", ggeFile);
    exit(1);
  }
  if (SDDS_CheckColumn(&SDDSin, "CnmS0", NULL, SDDS_ANY_FLOATING_TYPE, stderr) != SDDS_CHECK_OK) {
    fprintf(stderr, "Unable to find floating-point column \"CnmS0\" in file %s\n", ggeFile);
    exit(1);
  }
  if (SDDS_CheckColumn(&SDDSin, "CnmS2", NULL, SDDS_ANY_FLOATING_TYPE, stderr) != SDDS_CHECK_OK) {
    fprintf(stderr, "Unable to find floating-point column \"CnmS2\" in file %s\n", ggeFile);
    exit(1);
  }

  // Read first page
  if (SDDS_ReadPage(&SDDSin) <= 0) {
    fprintf(stderr, "Unable read the first page of %s\n", ggeFile);
    exit(1);
  }
  if ((rows = SDDS_RowCount(&SDDSin)) <= 1) {
    fprintf(stderr, "Too few z values in file %s\n", ggeFile);
    exit(1);
  }
  if (!(zValues = SDDS_GetColumnInDoubles(&SDDSin, "z"))) {
    fprintf(stderr, "Problem reading column z from %s\n", ggeFile);
    exit(1);
  }
  if (!(CnmS0Values = SDDS_GetColumnInDoubles(&SDDSin, "CnmS0"))) {
    fprintf(stderr, "Problem reading column CnmS0 from %s\n", ggeFile);
    exit(1);
  }
  if (!(CnmS2Values = SDDS_GetColumnInDoubles(&SDDSin, "CnmS2"))) {
    fprintf(stderr, "Problem reading column CnmS2 from %s\n", ggeFile);
    exit(1);
  }
  *Nz = rows;
  *z = zValues;
  *ggeD = CnmS0Values;
  *ggeS = CnmS2Values;
  for (i = 0; i < rows; i++) {
    (*ggeD)[i] *= -1;
    (*ggeS)[i] *= .25;
  }

  // Read second page
  if (SDDS_ReadPage(&SDDSin) <= 0) {
    fprintf(stderr, "Unable read the second page of %s\n", ggeFile);
    exit(1);
  }
  if ((rows = SDDS_RowCount(&SDDSin)) != *Nz) {
    fprintf(stderr, "Unequal row count between pages %s\n", ggeFile);
    exit(1);
  }
  if (!(zValues = SDDS_GetColumnInDoubles(&SDDSin, "z"))) {
    fprintf(stderr, "Problem reading column z from %s\n", ggeFile);
    exit(1);
  }
  if (!(CnmS0Values = SDDS_GetColumnInDoubles(&SDDSin, "CnmS0"))) {
    fprintf(stderr, "Problem reading column CnmS0 from %s\n", ggeFile);
    exit(1);
  }
  *ggeQ = CnmS0Values;
  for (i = 0; i < rows; i++) {
    if ((*z)[i] != zValues[i]) {
      fprintf(stderr, "z values differ between pages %s\n", ggeFile);
      exit(1);
    }
    (*ggeQ)[i] *= -2;
  }
  free(zValues);

  // Read third page
  if (SDDS_ReadPage(&SDDSin) <= 0) {
    fprintf(stderr, "Unable read the third page of %s\n", ggeFile);
    exit(1);
  }
  if ((rows = SDDS_RowCount(&SDDSin)) != *Nz) {
    fprintf(stderr, "Unequal row count between pages %s\n", ggeFile);
    exit(1);
  }
  if (!(zValues = SDDS_GetColumnInDoubles(&SDDSin, "z"))) {
    fprintf(stderr, "Problem reading column z from %s\n", ggeFile);
    exit(1);
  }
  if (!(CnmS0Values = SDDS_GetColumnInDoubles(&SDDSin, "CnmS0"))) {
    fprintf(stderr, "Problem reading column CnmS0 from %s\n", ggeFile);
    exit(1);
  }
  if (rows != *Nz) {
    fprintf(stderr, "Unequal row count between pages %s\n", ggeFile);
    exit(1);
  }
  for (i = 0; i < rows; i++) {
    if ((*z)[i] != zValues[i]) {
      fprintf(stderr, "z values differ between pages %s\n", ggeFile);
      exit(1);
    }
    (*ggeS)[i] -= 6.0 * CnmS0Values[i];
  }
  free(zValues);
  free(CnmS0Values);
  if (!SDDS_Terminate(&SDDSin)) {
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
    exit(1);
  }

  return;
}

void readInTrajectoryInput(char *input, double **x, double **xp, double **By, double *dz, int64_t *rows,
                           double *pCentral, int *Nsegments, double **zMagnetRef) {
  SDDS_DATASET SDDSin;
  int64_t i;
  int32_t segments;
  double *px, *pz, *z, dz0;
  char name[30];
  int ip;

  if (!SDDS_InitializeInput(&SDDSin, input)) {
    fprintf(stderr, "Unable to read file %s\n", input);
    exit(1);
  }
  if (SDDS_CheckParameter(&SDDSin, "pCentral", NULL, SDDS_ANY_FLOATING_TYPE, stderr) != SDDS_CHECK_OK) {
    fprintf(stderr, "Unable to find floating-point parameter \"pCentral\" in file %s\n", input);
    exit(1);
  }
  if (SDDS_CheckParameter(&SDDSin, "nMagnetSegments", NULL, SDDS_ANY_INTEGER_TYPE, stderr) != SDDS_CHECK_OK) {
    fprintf(stderr, "Unable to find integer paramenter \"nMagnetSegments\" in file %s\n", input);
    exit(1);
  }
  if (SDDS_CheckColumn(&SDDSin, "x", "m", SDDS_ANY_FLOATING_TYPE, stderr) != SDDS_CHECK_OK) {
    fprintf(stderr, "Unable to find floating-point column \"x\" with units \"m\" in file %s\n", input);
    exit(1);
  }
  if (SDDS_CheckColumn(&SDDSin, "z", "m", SDDS_ANY_FLOATING_TYPE, stderr) != SDDS_CHECK_OK) {
    fprintf(stderr, "Unable to find floating-point column \"x\" with units \"m\" in file %s\n", input);
    exit(1);
  }
  if (SDDS_CheckColumn(&SDDSin, "By", "T", SDDS_ANY_FLOATING_TYPE, stderr) != SDDS_CHECK_OK) {
    fprintf(stderr, "Unable to find floating-point column \"By\" with units \"T\" in file %s\n", input);
    exit(1);
  }
  if (SDDS_CheckColumn(&SDDSin, "px", NULL, SDDS_ANY_FLOATING_TYPE, stderr) != SDDS_CHECK_OK) {
    fprintf(stderr, "Unable to find floating-point column \"px\" in file %s\n", input);
    exit(1);
  }
  if (SDDS_CheckColumn(&SDDSin, "pz", NULL, SDDS_ANY_FLOATING_TYPE, stderr) != SDDS_CHECK_OK) {
    fprintf(stderr, "Unable to find floating-point column \"pz\" in file %s\n", input);
    exit(1);
  }
  if (SDDS_ReadPage(&SDDSin) <= 0) {
    fprintf(stderr, "Unable read the first page of %s\n", input);
    exit(1);
  }
  if ((*rows = SDDS_RowCount(&SDDSin)) <= 1) {
    fprintf(stderr, "Too few rows in file %s\n", input);
    exit(1);
  }
  if (!SDDS_GetParameterAsDouble(&SDDSin, "pCentral", pCentral)) {
    fprintf(stderr, "Problem reading parameter pCentral from %s\n", input);
    exit(1);
  }
  if (!SDDS_GetParameterAsLong(&SDDSin, "nMagnetSegments", &segments)) {
    fprintf(stderr, "Problem reading parameter nMagnetSegments from %s\n", input);
    exit(1);
  }
  *Nsegments = segments;
  for (ip = 1; ip <= *Nsegments; ip++) {
    sprintf(name, "zReference%d", ip);
    if (!SDDS_GetParameterAsDouble(&SDDSin, name, &(*zMagnetRef)[ip - 1])) {
      fprintf(stderr, "Problem reading parameter %s from %s\n", name, input);
      exit(1);
    }
  }
  if (!SDDS_GetParameterAsDouble(&SDDSin, "ZEntry", &zEntry) ||
      !SDDS_GetParameterAsDouble(&SDDSin, "ZVertex", &zVertex) ||
      !SDDS_GetParameterAsDouble(&SDDSin, "ZExit", &zExit)) {
    fprintf(stderr, "Problem reading parameters ZEntry, ZVertex, and ZExit from %s\n", input);
    exit(1);
  }
  if (!SDDS_GetParameterAsDouble(&SDDSin, "XEntry", &xEntry) ||
      !SDDS_GetParameterAsDouble(&SDDSin, "XVertex", &xVertex) ||
      !SDDS_GetParameterAsDouble(&SDDSin, "XExit", &xExit)) {
    fprintf(stderr, "Problem reading parameters XEntry, XVertex, and XExit from %s\n", input);
    exit(1);
  }

  if (!(*x = SDDS_GetColumnInDoubles(&SDDSin, "x"))) {
    fprintf(stderr, "Problem reading column x from %s\n", input);
    exit(1);
  }
  if (!(*By = SDDS_GetColumnInDoubles(&SDDSin, "By"))) {
    fprintf(stderr, "Problem reading column By from %s\n", input);
    exit(1);
  }
  if (!(z = SDDS_GetColumnInDoubles(&SDDSin, "z"))) {
    fprintf(stderr, "Problem reading column z from %s\n", input);
    exit(1);
  }
  if ((*dz = (z[*rows - 1] - z[0]) / (*rows - 1)) <= 0)
    bomb("dz is non-positive in trajectory file", NULL);
  for (i = 0; i < (*rows - 1); i++)
    if (fabs((dz0 = z[i + 1] - z[i]) / (*dz) - 1) > 1e-6)
      bomb("fractional dz variation of more than 1e-6 seen in trajectory file", NULL);
  if (!(px = SDDS_GetColumnInDoubles(&SDDSin, "px"))) {
    fprintf(stderr, "Problem reading column px from %s\n", input);
    exit(1);
  }
  if (!(pz = SDDS_GetColumnInDoubles(&SDDSin, "pz"))) {
    fprintf(stderr, "Problem reading column pz from %s\n", input);
    exit(1);
  }
  if (!(*xp = malloc(sizeof(**xp) * (*rows)))) {
    fprintf(stderr, "Problem allocating memory for xp\n");
  }
  for (i = 0; i < *rows; i++)
    (*xp)[i] = px[i] / pz[i];
  free(px);
  free(pz);
  if (!SDDS_Terminate(&SDDSin)) {
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
    exit(1);
  }
}

/*
  int refineNextMaximum(double *ggeD, int *zMaxInt, int edgeNum)
  {
  double max;
  int search, searchRet, count = 0;

  search = zMaxInt[edgeNum];
  max = ggeD[search];
  search++;
  do
  {
  if (fabs(ggeD[search]) > max)
  {
  max = fabs(ggeD[search]);
  searchRet = search;
  count = 0;
  }
  else
  count++; // prevent local mimima due to noise
  search++;
  }
  while (count < 10);

  return (searchRet);
  }
*/

/*** integrates integrand using the trapezoidal rule from startPt to endPt ***/
double integrateTrap(double *integrand, int startPt, int endPt, double dz) {
  double integral = 0.0;
  int ip;

  for (ip = startPt; ip < endPt; ip++)
    integral += 0.5 * (integrand[ip] + integrand[ip + 1]);
  integral *= dz;

  return (integral);
}

/*** sets heaviside step function such that the integrated bending field between  ***/
/*** zMaxInt[edgeNum] and zEdgeInt[edgeNum] matches that of field map = checkInt  ***/
void setStepFunction(double *heavside, double *maxD, int *zMaxInt, int *zEdgeInt, int edgeNum, double checkInt, double dz) {
  double integral;

  int ip, test = -1;

  while (test < 0) {
    /* set constant value in segment */
    for (ip = zMaxInt[edgeNum]; ip < zEdgeInt[edgeNum]; ip++)
      heavside[ip] = maxD[edgeNum];
    /* set endpoint to zero  */
    heavside[zEdgeInt[edgeNum]] = 0.0;
    /* set constant value in next segment */
    for (ip = zEdgeInt[edgeNum] + 1; ip <= zMaxInt[edgeNum + 1]; ip++)
      heavside[ip] = maxD[edgeNum + 1];
    integral = integrateTrap(heavside, zMaxInt[edgeNum], zMaxInt[edgeNum + 1], dz);
    /* adjust endpoint in present segment to match integral */
    heavside[zEdgeInt[edgeNum]] = (checkInt - integral) / dz;

    if (maxD[edgeNum] > maxD[edgeNum + 1]) { /* Assumes field in next segment is smaller */
      if ((maxD[edgeNum] > heavside[zEdgeInt[edgeNum]]) && (heavside[zEdgeInt[edgeNum]] > maxD[edgeNum + 1]))
        test = 1; /* endpoint correction is bewteen values in 2 segments -> success */
      else        /* Move edge location by one step and return to top */
        if (heavside[zEdgeInt[edgeNum]] > maxD[edgeNum])
          zEdgeInt[edgeNum]++;
        else
          zEdgeInt[edgeNum]--;
    } else { /* Assumes field in next segment is larger */
      if ((maxD[edgeNum + 1] > heavside[zEdgeInt[edgeNum]]) && (heavside[zEdgeInt[edgeNum]] > maxD[edgeNum]))
        test = 1; /* endpoint correction is bewteen values in 2 segments -> success */
      else        /* Move edge location by one step and return to top */
        if (heavside[zEdgeInt[edgeNum]] > maxD[edgeNum + 1])
          zEdgeInt[edgeNum]--;
        else
          zEdgeInt[edgeNum]++;
    }
  }
}

/*** Calculates dipole fringe field contributions ***/
FRINGE_INT3 computeDipoleFringeInt(double *ggeD, double *heavside, double *z, double *zEdge,
                                   int *zMaxInt, int edgeNum, double invRigidity, double dz) {
  FRINGE_INT3 dipoleFringeInt;
  double intK0 = 0.0;
  double intK1 = 0.0;
  double intK2 = 0.0;
  double sumHeav = heavside[zMaxInt[edgeNum]] + heavside[zMaxInt[edgeNum + 1]];
  double prodHeav = heavside[zMaxInt[edgeNum]] * heavside[zMaxInt[edgeNum + 1]];

  int ip;

  for (ip = zMaxInt[edgeNum]; ip < zMaxInt[edgeNum + 1]; ip++) {
    intK0 += 0.5 * ((z[ip] - zEdge[edgeNum]) * (heavside[ip] - ggeD[ip]) + (z[ip + 1] - zEdge[edgeNum]) * (heavside[ip + 1] - ggeD[ip + 1]));
    intK1 += 0.5 * ((ggeD[ip] - heavside[ip]) + (ggeD[ip + 1] - heavside[ip + 1]));
    intK2 += 0.5 * ((sumHeav - ggeD[ip]) * ggeD[ip] + (sumHeav - ggeD[ip + 1]) * ggeD[ip + 1] - 2.0 * prodHeav);
  }
  dipoleFringeInt.int1 = dz * invRigidity * intK1;
  dipoleFringeInt.int2 = dz * invRigidity * intK0;
  dipoleFringeInt.int3 = dz * invRigidity * invRigidity * intK2;

  return (dipoleFringeInt);
}

/*** Calculates quadrupole fringe field contributions ***/
FRINGE_INT3 computeQuadrupoleFringeInt(double *ggeQ, double *stepFuncQ, double *z, double *zEdge,
                                       int *zMaxInt, int edgeNum, double invRigidity, double dz) {
  FRINGE_INT3 quadFringeInt;
  double intQ0 = 0.0;
  double intQ1 = 0.0;
  double intQ2 = 0.0;

  int ip;

  for (ip = zMaxInt[edgeNum]; ip < zMaxInt[edgeNum + 1]; ip++) {
    intQ0 += 0.5 * ((ggeQ[ip] - stepFuncQ[ip]) + (ggeQ[ip + 1] - stepFuncQ[ip + 1]));
    intQ1 += 0.5 * ((z[ip] - zEdge[edgeNum]) * (ggeQ[ip] - stepFuncQ[ip]) + (z[ip + 1] - zEdge[edgeNum]) * (ggeQ[ip + 1] - stepFuncQ[ip + 1]));

    intQ2 += 0.5 * ((z[ip] - zEdge[edgeNum]) * (z[ip] - zEdge[edgeNum]) * (ggeQ[ip] - stepFuncQ[ip]) + (z[ip + 1] - zEdge[edgeNum]) * (z[ip + 1] - zEdge[edgeNum]) * (ggeQ[ip + 1] - stepFuncQ[ip + 1]));
  }
  quadFringeInt.int1 = dz * invRigidity * intQ0;
  quadFringeInt.int2 = dz * invRigidity * intQ1;
  quadFringeInt.int3 = dz * invRigidity * intQ2;

  return (quadFringeInt);
}

/*** Calculates sextupole/curvature ringe field contributions ***/
FRINGE_INT3 computeSextupoleFringeInt(double *ggeS, double *stepFuncS, double *z, double *zEdge,
                                      int *zMaxInt, int edgeNum, double invRigidity, double dz) {
  FRINGE_INT3 sextFringeInt;
  double intK4 = 0.0;
  double intK5 = 0.0;
  double intK6 = 0.0;

  int ip;

  for (ip = zMaxInt[edgeNum]; ip < zMaxInt[edgeNum + 1]; ip++) {
    intK4 += 0.5 * ((z[ip] - zEdge[edgeNum]) * (z[ip] - zEdge[edgeNum]) * (ggeS[ip] - stepFuncS[ip]) + (z[ip + 1] - zEdge[edgeNum]) * (z[ip + 1] - zEdge[edgeNum]) * (ggeS[ip + 1] - stepFuncS[ip + 1]));
    intK5 += 0.5 * ((z[ip] - zEdge[edgeNum]) * (ggeS[ip] - stepFuncS[ip]) + (z[ip + 1] - zEdge[edgeNum]) * (ggeS[ip + 1] - stepFuncS[ip + 1]));
    intK6 += 0.5 * ((ggeS[ip] - stepFuncS[ip]) + (ggeS[ip + 1] - stepFuncS[ip + 1]));
  }
  sextFringeInt.int1 = dz * invRigidity * intK6;
  sextFringeInt.int2 = dz * invRigidity * intK5;
  sextFringeInt.int3 = dz * invRigidity * intK4;

  return (sextFringeInt);
}

FRINGE_INT3 computeDipSextFringeInt(double *ggeD, double *stepFuncD, double *ggeS, double *stepFuncS,
                                    double *z, double *zEdge, int *zMaxInt, int edgeNum, double invRigidity, double dz) {
  FRINGE_INT3 dipSextFringeInt;

  double intI5 = 0.0;
  double temp6a, temp6b, intI6 = 0.0;
  double temp7, intI7 = 0.0;

  int ip, id;

  for (ip = zMaxInt[edgeNum]; ip < zMaxInt[edgeNum + 1]; ip++) {
    intI5 += 0.5 * invRigidity * (z[ip] - zEdge[edgeNum]) * (z[ip] - zEdge[edgeNum]) * (ggeS[ip] - stepFuncS[ip]) * ggeD[ip];

    temp6a = 0.5 * dz * invRigidity * (ggeD[ip] - stepFuncD[ip]);
    temp6b = 0.5 * dz * invRigidity * (ggeS[ip] - stepFuncS[ip]);

    temp7 = 0.5 * dz * invRigidity * (ggeS[ip] - stepFuncS[ip]);
    for (id = zMaxInt[edgeNum]; id < ip - 1; id++) {
      intI6 += 0.5 * (z[ip] - z[id]) * (temp6a * (ggeS[id] - stepFuncS[id]) + temp6b * (ggeD[id] - stepFuncD[id]));
      intI6 += 0.5 * (z[ip] - z[id + 1]) * (temp6a * (ggeS[id + 1] - stepFuncS[id + 1]) + temp6b * (ggeD[id + 1] - stepFuncD[id + 1]));

      intI7 += 0.5 * temp7 * (z[ip] - z[id]) * ggeD[id];
      intI7 += 0.5 * temp7 * (z[ip] - z[id + 1]) * ggeD[id + 1];
    }
    intI5 += 0.5 * invRigidity * (z[ip + 1] - zEdge[edgeNum]) * (z[ip + 1] - zEdge[edgeNum]) * (ggeS[ip + 1] - stepFuncS[ip + 1]) * ggeD[ip + 1];

    temp6a = 0.5 * dz * invRigidity * (ggeD[ip + 1] - stepFuncD[ip + 1]);
    temp6b = 0.5 * dz * invRigidity * (ggeS[ip + 1] - stepFuncS[ip + 1]);

    temp7 = 0.5 * dz * invRigidity * (ggeS[ip + 1] - stepFuncS[ip + 1]);
    for (id = zMaxInt[edgeNum]; id < ip; id++) {
      intI6 += 0.5 * (z[ip + 1] - z[id]) * (temp6a * (ggeS[id] - stepFuncS[id]) + temp6b * (ggeD[id] - stepFuncD[id]));
      intI6 += 0.5 * (z[ip + 1] - z[id + 1]) * (temp6a * (ggeS[id + 1] - stepFuncS[id + 1]) + temp6b * (ggeD[id + 1] - stepFuncD[id + 1]));

      intI7 += 0.5 * temp7 * (z[ip + 1] - z[id]) * ggeD[id];
      intI7 += 0.5 * temp7 * (z[ip + 1] - z[id + 1]) * ggeD[id + 1];
    }
  }
  dipSextFringeInt.int1 = 0.5 * dz * invRigidity * intI5;
  dipSextFringeInt.int2 = 0.5 * dz * invRigidity * intI6;
  dipSextFringeInt.int3 = dz * invRigidity * intI7;

  return (dipSextFringeInt);
}

void readInAndMoveGGE(char *ggeFile, double **z, double **ggeD, double **ggeQ, double **ggeS, int *Nz, double xEntry) {
  SDDS_TABLE SDDS_inputGGE;

  GGE_m_ell **gges, **deriv_gges;
  double *zValues;

  int *m;

  char ggeName[20];
  int readCode;
  int il, l_2Max = 0;
  int im, mMax = 0;
  int iz, rows = 0, lastRows = 0;

  if (!SDDS_InitializeInput(&SDDS_inputGGE, ggeFile)) {
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
    exit(1);
  }

  l_2Max = 0;
  sprintf(ggeName, "CnmS%i", 2 * l_2Max);
  while (SDDS_GetColumnIndex(&SDDS_inputGGE, ggeName) >= 0) {
    l_2Max++;
    sprintf(ggeName, "CnmS%i", 2 * l_2Max);
  }
  l_2Max -= 1;

  // find number of gradients
  mMax = 0;
  while ((readCode = SDDS_ReadTable(&SDDS_inputGGE)) > 0) {
    mMax++;
    if ((rows = SDDS_RowCount(&SDDS_inputGGE)) <= 1)
      SDDS_Bomb("Too few points in z in GGE file");
    if (readCode > 1 && lastRows != rows)
      SDDS_Bomb("Inconsistent number of z points in GGE file");
    lastRows = rows;
  }

  if (!SDDS_Terminate(&SDDS_inputGGE)) {
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
    exit(1);
  }

  // allocate memory
  gges = calloc(mMax, sizeof(GGE_m_ell *));
  deriv_gges = calloc(mMax, sizeof(GGE_m_ell *));
  for (im = 0; im < mMax; im++) {
    gges[im] = calloc(l_2Max + 1, sizeof(GGE_m_ell));
    deriv_gges[im] = calloc(l_2Max + 1, sizeof(GGE_m_ell));
    for (il = 0; il <= l_2Max; il++) {
      gges[im][il].grad = calloc(rows, sizeof(double));
      deriv_gges[im][il].grad = calloc(rows, sizeof(double));
    }
  }
  m = calloc(mMax, sizeof(int));

  if (!SDDS_InitializeInput(&SDDS_inputGGE, ggeFile)) {
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
    exit(1);
  }

  // Initialize gradients
  im = 0;
  zValues = NULL; // suppress compiler warning
  while ((readCode = SDDS_ReadTable(&SDDS_inputGGE)) > 0) {
    SDDS_GetParameter(&SDDS_inputGGE, "m", &m[im]);

    if (!(zValues = SDDS_GetColumnInDoubles(&SDDS_inputGGE, "z"))) {
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
      exit(1);
    }

    for (il = 0; il <= l_2Max; il++) {
      snprintf(ggeName, sizeof(ggeName), "CnmS%i", 2 * il);
      if (!(gges[im][il].grad = SDDS_GetColumnInDoubles(&SDDS_inputGGE, ggeName))) {
        SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
        exit(1);
      }
    }
    for (il = 0; il <= l_2Max; il++) {
      snprintf(ggeName, sizeof(ggeName), "dCnmS%i/dz", 2 * il);
      if (!(deriv_gges[im][il].grad = SDDS_GetColumnInDoubles(&SDDS_inputGGE, ggeName))) {
        SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
        exit(1);
      }
    }
    im++;
  }

  // Use the generalized gradients to compute Bnormal on a cylinder centered
  // at x=xEntry with radius=xEntry (could be generalized)
  double radius = xEntry;
  int iangle, Nangle = 256;
  double **Bnorm;
  double *dipole, *DdipoleD2, *quadrupole, *sextupole;
  Bnorm = calloc(Nangle, sizeof(double *));
  for (iangle = 0; iangle < Nangle; iangle++)
    Bnorm[iangle] = calloc(rows, sizeof(double));
  computeBOnCylinder(xEntry, radius, Bnorm, m, gges, mMax, l_2Max, rows, Nangle);
  dipole = calloc(rows, sizeof(double));
  DdipoleD2 = calloc(rows, sizeof(double));
  quadrupole = calloc(rows, sizeof(double));
  sextupole = calloc(rows, sizeof(double));

  // Use Bnormal on new cylinder to compute GGEs needed for calculation along new axis
  computeRequiredGGEs(dipole, DdipoleD2, quadrupole, sextupole, Bnorm, radius, zValues[1] - zValues[0], Nangle, rows);
  for (iangle = 0; iangle < Nangle; iangle++)
    free(Bnorm[iangle]);
  free(Bnorm);

  // Set gradients etc.
  *Nz = rows;
  *z = zValues;
  *ggeD = dipole;
  *ggeQ = quadrupole;
  *ggeS = sextupole;
  for (iz = 0; iz < rows; iz++) {
    (*ggeD)[iz] *= -1.0;
    (*ggeQ)[iz] *= -2.0;
    (*ggeS)[iz] *= -6.0;
    (*ggeS)[iz] += 0.25 * DdipoleD2[iz];
  }

  for (im = 0; im < mMax; im++) {
    free(gges[im]);
    free(deriv_gges[im]);
  }
  free(gges);
  free(deriv_gges);
  free(m);

  /*
    These arrays are used to return the results to the caller
    free(dipole);
    free(DdipoleD2);
    free(quadrupole);
    free(sextupole);
  */

  return;
}
void computeBOnCylinder(double xCenter, double radius, double **Bnorm, int *m,
                        GGE_m_ell **gges, int mMax, int l_2Max, int Nz, int Nangle) {
  double r, x, y, theta, phi;
  double Br, Bphi, Bx, By;
  double dAngle = 2.0 * PI / (double)Nangle;

  int iz, iangle, im, il;

  for (iz = 0; iz < Nz; iz++)
    for (iangle = 0; iangle < Nangle; iangle++) {
      theta = dAngle * (double)iangle;
      x = radius * cos(theta) + xCenter;
      y = radius * sin(theta);
      r = sqrt(x * x + y * y);
      phi = atan2(y, x);

      Br = Bphi = 0.0;
      for (im = 0; im < mMax; im++) {
        double mfact, term, sin_mphi, cos_mphi;
        mfact = (double)factorial(m[im]);
        sin_mphi = sin(m[im] * phi);
        cos_mphi = cos(m[im] * phi);

        for (il = 0; il <= l_2Max; il++) {
          term = gges[im][il].grad[iz] * intPower(-1.0, il) * mfact / (intPower(2.0, 2 * il) * factorial(il) * factorial(il + m[im])) * intPower(r, 2 * il + m[im] - 1);

          Br += term * (2 * il + m[im]) * sin_mphi;
          Bphi += m[im] * term * cos_mphi;
        }
      }
      Bx = (Br * cos(phi) - Bphi * sin(phi));
      By = (Br * sin(phi) + Bphi * cos(phi));

      Bnorm[iangle][iz] = Bx * cos(theta) + By * sin(theta);
    }
  return;
}

void computeRequiredGGEs(double *dipole, double *DdipoleD2, double *quadrupole, double *sextupole,
                         double **Bnorm, double rho, double dz, int Nangle, int Nz) {
  COMPLEX **BnormFFT;
  COMPLEX *Btemp, *gge;

  double *k;

  double dk, factor;

  int iz, ik, iangle, im, il;

  BnormFFT = calloc(Nangle, sizeof(COMPLEX *));
  for (iangle = 0; iangle < Nangle; iangle++)
    BnormFFT[iangle] = calloc(Nz, sizeof(COMPLEX));
  Btemp = calloc(Nangle, sizeof(COMPLEX));
  /* Take FFT of Brho values vs phi for each z plane */
  for (iz = 0; iz < Nz; iz++) {
    for (iangle = 0; iangle < Nangle; iangle++) {
      Btemp[iangle].re = Bnorm[iangle][iz];
      Btemp[iangle].im = 0.0;
    }
    FFT(Btemp, -1, Nangle);
    for (iangle = 0; iangle < Nangle; iangle++) {
      BnormFFT[iangle][iz].re = Btemp[iangle].re;
      BnormFFT[iangle][iz].im = Btemp[iangle].im;
    }
  }
  free(Btemp);

  /* Take FFT of Brho values vs z for lowest harmonics */
  for (iangle = 1; iangle <= 3; iangle++)
    FFT(BnormFFT[iangle], -1, Nz);

  dk = 2.0 * PI / (dz * (double)Nz);
  k = calloc(Nz, sizeof(double));
  for (ik = 0; ik < Nz / 2; ik++)
    k[ik] = dk * ik;
  k[0] = 1e-12 * dk;
  for (ik = Nz / 2; ik < Nz; ik++)
    k[ik] = -dk * (Nz - ik);

  gge = calloc(Nz, sizeof(COMPLEX));

  // Use the formulas to compute C_{10}, C_{12}, C_{20}, and C_{30}
  il = 0;
  im = 1;
  for (ik = 0; ik < Nz; ik++) {
    factor = ipow(-1, il / 2) * ipow(k[ik], im + il - 1) / Imp(k[ik] * rho, im) / (ipow(2, im) * dfactorial(im)) / (0.5 * Nz * Nangle);
    gge[ik].re = -BnormFFT[im][ik].im * factor;
    gge[ik].im = BnormFFT[im][ik].re * factor;
  }
  FFT(gge, 1, Nz);
  for (ik = 0; ik < Nz; ik++)
    dipole[ik] = gge[ik].re;
  il = 2;
  for (ik = 0; ik < Nz; ik++) {
    factor = ipow(-1, il / 2) * ipow(k[ik], im + il - 1) / Imp(k[ik] * rho, im) / (ipow(2, im) * dfactorial(im)) / (0.5 * Nz * Nangle);
    gge[ik].re = -BnormFFT[im][ik].im * factor;
    gge[ik].im = BnormFFT[im][ik].re * factor;
  }
  FFT(gge, 1, Nz);
  for (ik = 0; ik < Nz; ik++)
    DdipoleD2[ik] = gge[ik].re;

  im = 2;
  il = 0;
  for (ik = 0; ik < Nz; ik++) {
    factor = ipow(-1, il / 2) * ipow(k[ik], im + il - 1) / Imp(k[ik] * rho, im) / (ipow(2, im) * dfactorial(im)) / (0.5 * Nz * Nangle);
    gge[ik].re = -BnormFFT[im][ik].im * factor;
    gge[ik].im = BnormFFT[im][ik].re * factor;
  }
  FFT(gge, 1, Nz);
  for (ik = 0; ik < Nz; ik++)
    quadrupole[ik] = gge[ik].re;

  im = 3;
  for (ik = 0; ik < Nz; ik++) {
    factor = ipow(-1, il / 2) * ipow(k[ik], im + il - 1) / Imp(k[ik] * rho, im) / (ipow(2, im) * dfactorial(im)) / (0.5 * Nz * Nangle);
    gge[ik].re = -BnormFFT[im][ik].im * factor;
    gge[ik].im = BnormFFT[im][ik].re * factor;
  }
  FFT(gge, 1, Nz);
  for (ik = 0; ik < Nz; ik++)
    sextupole[ik] = gge[ik].re;

  for (iangle = 0; iangle < Nangle; iangle++)
    free(BnormFFT[iangle]);
  free(BnormFFT);
  free(gge);
  free(k);

  return;
}

/*******************************************/
/* Calculates x^n assuming n is an integer */
/*******************************************/
double intPower(double x, int n) {
  double out = 1.0;
  int j;

  for (j = 0; j < abs(n); j++)
    out = out * x;

  if (n < 0)
    out = 1.0 / out;

  return (out);
}

void FFT(COMPLEX *field, int32_t isign, int32_t npts) {
  double *real_imag;
  int32_t i;

  real_imag = tmalloc(sizeof(double) * (2 * npts + 2));
  for (i = 0; i < npts; i++) {
    real_imag[2 * i] = field[i].re;
    real_imag[2 * i + 1] = field[i].im;
  }
  if (isign == -1) {
    complexFFT(real_imag, npts, 0);
    for (i = 0; i < npts; i++) {
      field[i].re = npts * real_imag[2 * i];
      field[i].im = npts * real_imag[2 * i + 1];
    }
  } else {
    complexFFT(real_imag, npts, INVERSE_FFT);
    for (i = 0; i < npts; i++) {
      field[i].re = real_imag[2 * i];
      field[i].im = real_imag[2 * i + 1];
    }
  }
  free(real_imag);
}

double Imp(double kR, long m) {
  double t1;

  // For some reason GSL didn't work for all, so I wrote my own code to compute I_n(x)
  //  t1 = (BesIn(fabs(kR), m-1) + BesIn(fabs(kR), m+1))/2;
  t1 = (BesselFuncIn(fabs(kR), m - 1) + BesselFuncIn(fabs(kR), m + 1)) / 2;
  if (kR < 0)
    return t1 * ipow(-1, m + 1);
  return t1;
}

#ifdef DONT_USE_THIS
double BesIn(double x, long order) {
  double result;
  if (order == 0) {
    result = dbesi0(x);
  } else if (order == 1) {
    result = dbesi1(x);
  } else {
    if (order < 0)
      order = -order;
    /* Need this code to compensate for issues with gsl_sf_bessel_Inu() domain */
    if (x > 0)
      result = gsl_sf_bessel_In(order, x);
    else {
      if (((long)fabs(order)) % 2)
        result = -gsl_sf_bessel_In(order, -x);
      else
        result = gsl_sf_bessel_In(order, -x);
    }
  }
  return result;
}
#endif

/* Computes I_order(x) assuming x>0, but order can be any real number */
double BesselFuncIn(double x, double order) {
  double jDub, f = 1.0;
  double nu = order;

  int j;

  if (x < 8.0 + nu * nu / 12.0) {
    if (fabs(nu + (int)fabs(nu)) < 1.e-12)
      nu = -nu;
    j = (int)(x + 3.0);
    if ((nu - 5.0) * (nu + 10.0 + 1.5 * x) < 0.0)
      j = j + (int)(4.0 - 2.0 * nu / 3.0);
    do {
      f = 1.0 + f * x * x / (4.0 * (double)j * ((double)j + nu));
      j--;
    } while (j > 0);
    if (fabs(nu) < 1.e-12) {
      return (f);
    }
    f = f / nu;
    do {
      f = 2.0 * f * nu / x;
      nu += 1.0;
    } while (nu <= 3.0);
    jDub = 2.0 / (3.0 * nu * nu) - 1.0;
    jDub = 1.0 + 2.0 * jDub / (7.0 * nu * nu);
    jDub = jDub / (30.0 * nu * nu);
    jDub = (jDub - 1.0) / (12.0 * nu);
    jDub += nu * (1.0 - log(2.0 * nu / x));
    f = f * exp(jDub) * sqrt(nu / (2.0 * PI));
    return (f);
  } else {
    j = (int)(5.0 + 640.0 / (x * x) + 0.7 * fabs(nu));
    if (x >= 100.0)
      j = j - (int)((0.5 - sqrt(247.0 / (x * sqrt(x)))) * fabs(nu));
    do {
      jDub = ((double)j - 0.5) * ((double)j - 0.5);
      f = 1.0 + f * (jDub - nu * nu) / (2.0 * x * (double)j);
      j--;
    } while (j > 0);
    f = f * exp(x) / sqrt(2.0 * PI * x);
  }

  return (f);
}
