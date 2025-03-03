#include "mdb.h"
#include "SDDS.h"
#include "scan.h"
#include "fftpackC.h"
#include "gsl/gsl_sf_bessel.h"

#define TWOPI 6.28318530717958647692528676656

#if defined(linux) || (defined(_WIN32) && !defined(_MINGW))
#  include <omp.h>
#else
#  define NOTHREADS 1
#endif

typedef struct COMPLEX {
  double re;
  double im;
} COMPLEX;

typedef struct {
  double *Brho, *Bz; /* won't be changed once read */
  long Nphi, Nz;
  double dphi, dz, rho;
  double zMin;
} FIELDS_ON_BOUNDARY;

typedef struct {
  double *x, *y, *z;
  double *Bx, *By, *Bz;
  long n;
} FIELD_MAP;

typedef struct {
  double rms, mad, max;
  double fracRms, fracMad, fracMax;
} ALL_RESIDUALS;

typedef struct {
  short haveData;
  short skew;              /* if non-zero, these are skew terms */
  long nz;                 /* number of z points */
  double dz;               /* z spacing */
  double xCenter, yCenter; /* center of the expansion in magnet coordinate system */
  double xMax, yMax;       /* half-aperture of the field expansion in expansion coordinate system */
  double zMin, zMax;       /* minimum and maximum z values */
  long nm;                 /* number of values of m (angular harmonic) */
  long *m;                 /* value of m */
  long nGradients;         /* number of gradient functions per m */
  double ***Cmn;           /* generalized gradient: Cnms[im][in][iz] */
  double ***dCmn_dz;       /* z derivative of generalized gradient */
} BGGEXP_DATA;

int ReadInputFile(FIELDS_ON_BOUNDARY *fieldsOnBoundary, char *inputFile,
                  char *zName, char *phiName, char *BrhoName, char *BzName, char *rhoName);
void freeFieldsOnBoundary(FIELDS_ON_BOUNDARY *fob);
double BesIn(double x, long order);
int SetUpOutputFile(SDDS_DATASET *SDDSout, char *filename, long skew, long derivatives);
int StartPageLocal(SDDS_DATASET *SDDSout, FIELDS_ON_BOUNDARY *fob, long m, char *type);
void FFT(COMPLEX *field, int32_t isign, int32_t npts);

void readBGGExpData(BGGEXP_DATA *bggexpData, char *filename, char *nameFragment, short skew);
void freeBGGExpData(BGGEXP_DATA *bggexpData);
double evaluateGGEFit(FIELDS_ON_BOUNDARY *fob, BGGEXP_DATA bggexpData[2],
                      long derivatives, long multipoles, double significance, unsigned long flags,
                      long varyDerivatives, ALL_RESIDUALS *allResiduals);
void readFieldMap(char *fieldMapFile, FIELD_MAP *fmData);
int evaluateGGEAndOutput(char *outputFile, long nrho, long nphi, char *normalFile, char *skewFile,
                         FIELDS_ON_BOUNDARY *fob);
int computeGGE(FIELDS_ON_BOUNDARY *fieldsOnBoundary, char *normalOutput, char *skewOutput,
               long derivatives, long multipoles, long fundamental, long varyDerivatives);

#define SET_INPUTFILE 0
#define SET_NORMAL 1
#define SET_SKEW 2
#define SET_DERIVATIVES 3
#define SET_MULTIPOLES 4
#define SET_FUNDAMENTAL 5
#define SET_AUTO_TUNE 6
#define SET_EVALUATE 7
#define SET_THREADS 8
#define SET_VARY_DERIVATIVES 9
#define N_OPTIONS 10

char *option[N_OPTIONS] = {
  "input", "normal", "skew", "derivatives", "multipoles", "fundamental", "autotune", "evaluate", "threads",
  "varyderivatives"};

#define USAGE "computeCBGGE -input=<filename>[,z=<columnName>][,phi=<columnName>][,Brho=<columnName>][,rho=<parameterName>][,Bz=<columnName>]\n\
              -normal=<output> [-skew=<output>] [-threads=<integer>]\n\
              [-derivatives=<integer>] [-multipoles=<integer>] [-fundamental=<integer>]\n\
              [-varyDerivatives] [-evaluate=<filename>[,nrho=<integer>][,nphi=<integer>]\n \
              [-autotune=[,significance=<fieldValue>][,minimize={rms|mav|maximum}][,increaseOnly][,verbose][,log=<filename>]]\n\
-input       Single-page file giving (z, phi, Brho) on circular cylinder of radius rho.\n\
             If solenoidal fields are desired, file can also include Bz, but column must be named with Bz qualifier.\n\
-normal      Output file for normal-component generalized gradients.\n\
-skew        Output file for skew-component generalized gradients. If the input data\n\
             has non-zero Bz on axis, this option is essential.\n\
-threads     Number of threads to use for computations (default: 1).\n\
-derivatives Number of derivatives vs z desired in output. Default: 7\n\
-multipoles  Number of multipoles desired in output. Default: 8\n\
-fundamental Fundamental multipole of sequence. 0=none (default), 1=dipole, 2=quadrupole, etc.\n\
-varyDerivatives\n\
             If given, the number of derivatives used varies with multipole order so as to\n\
             maintain an approximately consistent maximum transverse order.\n\
             Recommended to try if GGE does not fit data well.\n\
-evaluate    Evaluate the GGE over the interior region, including the boundary, using the specified\n\
             number of radial and angular samples. Defaults to the nrho=1 and same angles as in input.\n\
-autotune    Seeks to minimize the number of multipoles and derivatives to avoid using terms\n\
             that do not contribute to a good fit at the given level of significance. The user can\n\
             choose to minimize the maximum error (default), the rms error, or the mean absolute value\n\
             error. The fit is evaluated with respect to the input file.\n\n\
Circular-cylinder Boundary Generalized Gradient Expansion by Michael Borland, Robert Soliday, and Ryan Lindberg."

#define AUTOTUNE_VERBOSE 0x0001UL
#define AUTOTUNE_RMS 0x0002UL
#define AUTOTUNE_MAXIMUM 0x0004UL
#define AUTOTUNE_MAV 0x0008UL
#define AUTOTUNE_EVALONLY 0x0010UL
#define AUTOTUNE_MODE_SET 0x0100UL
#define AUTOTUNE_ACTIVE 0x0200UL
#define AUTOTUNE_LOG 0x0400UL
#define AUTOTUNE_INCRONLY 0x0800UL
char *modeOption[3] = {"rms", "maximum", "mav"};

int threads = 1;

int main(int argc, char **argv) {
  SCANNED_ARG *scanned;
  long i_arg;
  long multipoles, derivatives, fundamental = 0;
  long maxMultipoles = 8, maxDerivatives = 7;
  char *inputFile = NULL, *normalOutputFile = NULL, *skewOutputFile = NULL;
  char *evaluationOutput = NULL;
  long evaluation_nRho = -1, evaluation_nPhi = -1;
  double autoTuneSignificance = 1e-12;
  FIELDS_ON_BOUNDARY fieldsOnBoundary;
  double bestResidual;
  long bestMultipoles, bestDerivatives;
  unsigned long autoTuneFlags = 0, dummyFlags;
  char *autoTuneModeString;
  char *BrhoName = "Brho", *zName = "z", *phiName = "phi", *rhoName = "rho", *BzName = NULL;
  BGGEXP_DATA bggexpData[2];
  SDDS_DATASET SDDS_autoTuneLog;
  char *autoTuneLogFile = NULL;
  long iAutoTuneLog = 0, varyDerivatives = 0;
  ALL_RESIDUALS allResiduals;
  long minMultipoles;

  /* Using this routine can prevent exceptions when large harmonics are used, but the data thus
     obtained is suspect.
     gsl_set_error_handler_off();
  */

  argc = scanargs(&scanned, argc, argv);
  if (argc < 2 || argc > (2 + N_OPTIONS)) {
    fprintf(stderr, "%s\n", USAGE);
    return (1);
  }
  for (i_arg = 1; i_arg < argc; i_arg++) {
    if (scanned[i_arg].arg_type == OPTION) {
      /* process options here */
      switch (match_string(scanned[i_arg].list[0], option, N_OPTIONS, 0)) {
      case SET_INPUTFILE:
        if (scanned[i_arg].n_items < 2) {
          fprintf(stderr, "invalid -input syntax\n%s\n", USAGE);
          return (1);
        }
        inputFile = scanned[i_arg].list[1];
        scanned[i_arg].n_items -= 2;
        if (!scanItemList(&dummyFlags, scanned[i_arg].list + 2, &scanned[i_arg].n_items, 0,
                          "z", SDDS_STRING, &zName, 1, 0,
                          "Brho", SDDS_STRING, &BrhoName, 1, 0,
                          "phi", SDDS_STRING, &phiName, 1, 0,
                          "rho", SDDS_STRING, &rhoName, 1, 0,
                          "Bz", SDDS_STRING, &BzName, 1, 0,
                          NULL)) {
          fprintf(stderr, "invalid -input syntax\n%s\n", USAGE);
          return 1;
        }
        break;
      case SET_NORMAL:
        if (scanned[i_arg].n_items != 2) {
          fprintf(stderr, "invalid -normal syntax\n%s\n", USAGE);
          return (1);
        }
        normalOutputFile = scanned[i_arg].list[1];
        break;
      case SET_SKEW:
        if (scanned[i_arg].n_items != 2) {
          fprintf(stderr, "invalid -skew syntax\n%s\n", USAGE);
          return (1);
        }
        skewOutputFile = scanned[i_arg].list[1];
        break;
      case SET_EVALUATE:
        if (scanned[i_arg].n_items < 2) {
          fprintf(stderr, "invalid -evaluate syntax\n%s\n", USAGE);
          return (1);
        }
        evaluationOutput = scanned[i_arg].list[1];
        scanned[i_arg].n_items -= 2;
        evaluation_nRho = -1;
        evaluation_nPhi = -1;
        if (scanned[i_arg].n_items > 0 &&
            !scanItemList(&dummyFlags, scanned[i_arg].list + 2, &scanned[i_arg].n_items, 0,
                          "nrho", SDDS_LONG, &evaluation_nRho, 1, 0,
                          "nphi", SDDS_LONG, &evaluation_nPhi, 1, 0,
                          NULL)) {
          fprintf(stderr, "invalid -evaluation syntax\n%s\n", USAGE);
          exit(1);
        }
        break;
      case SET_DERIVATIVES:
        if (scanned[i_arg].n_items != 2 ||
            sscanf(scanned[i_arg].list[1], "%ld", &maxDerivatives) != 1 ||
            maxDerivatives <= 0) {
          fprintf(stderr, "invalid -derivatives syntax\n%s\n", USAGE);
          return (1);
        }
        break;
      case SET_MULTIPOLES:
        if (scanned[i_arg].n_items != 2 ||
            sscanf(scanned[i_arg].list[1], "%ld", &maxMultipoles) != 1 ||
            maxMultipoles <= 0) {
          fprintf(stderr, "invalid -multipoles syntax\n%s\n", USAGE);
          return (1);
        }
        break;
      case SET_FUNDAMENTAL:
        if (scanned[i_arg].n_items != 2 ||
            sscanf(scanned[i_arg].list[1], "%ld", &fundamental) != 1 ||
            fundamental < 0) {
          fprintf(stderr, "invalid -fundamental syntax\n%s\n", USAGE);
          return (1);
        }
        break;
      case SET_AUTO_TUNE:
        autoTuneSignificance = 1e-12;
        autoTuneFlags = 0;
        scanned[i_arg].n_items -= 1;
        if (scanned[i_arg].n_items > 0 &&
            (!scanItemList(&autoTuneFlags, scanned[i_arg].list + 1, &scanned[i_arg].n_items, 0,
                           "verbose", -1, NULL, 0, AUTOTUNE_VERBOSE,
                           "increaseonly", -1, NULL, 0, AUTOTUNE_INCRONLY,
                           "evaluate", -1, NULL, 0, AUTOTUNE_EVALONLY,
                           "significance", SDDS_DOUBLE, &autoTuneSignificance, 1, 0,
                           "minimize", SDDS_STRING, &autoTuneModeString, 1, AUTOTUNE_MODE_SET,
                           "log", SDDS_STRING, &autoTuneLogFile, 1, AUTOTUNE_LOG,
                           NULL) ||
             autoTuneSignificance <= 0)) {
          fprintf(stderr, "invalid -autotune syntax\n%s\n", USAGE);
          return (1);
        }
        autoTuneFlags |= AUTOTUNE_ACTIVE;
        if (autoTuneFlags & AUTOTUNE_MODE_SET) {
          switch (match_string(autoTuneModeString, modeOption, 3, 0)) {
          case 0:
            autoTuneFlags |= AUTOTUNE_RMS;
            break;
          case 1:
            autoTuneFlags |= AUTOTUNE_MAXIMUM;
            break;
          case 2:
            autoTuneFlags |= AUTOTUNE_MAV;
            break;
          default:
            SDDS_Bomb("invalid mode for autotune minimization. Use rms or maximum.");
            break;
          }
        } else
          autoTuneFlags |= AUTOTUNE_MAXIMUM;
        break;
      case SET_THREADS:
        if (scanned[i_arg].n_items != 2 || sscanf(scanned[i_arg].list[1], "%d", &threads) != 1 || threads <= 0)
          SDDS_Bomb("invalid -threads syntax: give an value greater than 0");
        break;
      case SET_VARY_DERIVATIVES:
        varyDerivatives = 1;
        break;
      default:
        fprintf(stderr, "unknown option given\n%s\n", USAGE);
        return (1);
        break;
      }
    } else {
      fprintf(stderr, "Unreconized option\n%s\n", USAGE);
      return (1);
    }
  }

  if (inputFile == NULL) {
    fprintf(stderr, "%s\n", USAGE);
    return (1);
  }

#if !defined(NOTHREADS)
  omp_set_num_threads(threads);
#endif

  if (autoTuneFlags & AUTOTUNE_LOG) {
    if (SDDS_InitializeOutput(&SDDS_autoTuneLog, SDDS_BINARY, 1, NULL, "computeRBGGE autotune output", autoTuneLogFile) != 1 ||
        SDDS_DefineParameter(&SDDS_autoTuneLog, "OptimalMultipoles", "m$bopt$n", NULL, "Optimal number of multipoles",
                             NULL, SDDS_LONG, NULL) == -1 ||
        SDDS_DefineParameter(&SDDS_autoTuneLog, "OptimalDerivatives", "d$bopt$n", NULL, "Optimal number of derivatives",
                             NULL, SDDS_LONG, NULL) == -1 ||
        SDDS_DefineParameter(&SDDS_autoTuneLog, "OptimalResidual", "r$bopt$n", "T", "Optimal residual",
                             NULL, SDDS_DOUBLE, NULL) == -1 ||
        !SDDS_DefineSimpleParameter(&SDDS_autoTuneLog, "OptimumLabel", NULL, SDDS_STRING) ||
        SDDS_DefineColumn(&SDDS_autoTuneLog, "m", NULL, NULL, "Number of multipoles", NULL, SDDS_LONG, 0) == -1 ||
        SDDS_DefineColumn(&SDDS_autoTuneLog, "d", NULL, NULL, "Number of derivatives", NULL, SDDS_LONG, 0) == -1 ||
        SDDS_DefineSimpleColumn(&SDDS_autoTuneLog, "RmsError", "T", SDDS_DOUBLE) != 1 ||
        SDDS_DefineSimpleColumn(&SDDS_autoTuneLog, "MaximumError", "T", SDDS_DOUBLE) != 1 ||
        SDDS_DefineSimpleColumn(&SDDS_autoTuneLog, "MadError", "T", SDDS_DOUBLE) != 1 ||
        SDDS_DefineSimpleColumn(&SDDS_autoTuneLog, "FractionalRmsError", NULL, SDDS_DOUBLE) != 1 ||
        SDDS_DefineSimpleColumn(&SDDS_autoTuneLog, "FractionalMaximumError", NULL, SDDS_DOUBLE) != 1 ||
        SDDS_DefineSimpleColumn(&SDDS_autoTuneLog, "FractionalMadError", NULL, SDDS_DOUBLE) != 1 ||
        !SDDS_WriteLayout(&SDDS_autoTuneLog) ||
        !SDDS_StartPage(&SDDS_autoTuneLog, maxMultipoles * maxDerivatives)) {
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
      return (1);
    }
  }

  bestResidual = DBL_MAX;
  bestDerivatives = maxDerivatives;
  bestMultipoles = maxMultipoles;

  memset(&fieldsOnBoundary, 0, sizeof(fieldsOnBoundary));
  if (ReadInputFile(&fieldsOnBoundary, inputFile, zName, phiName, BrhoName, BzName, rhoName))
    SDDS_Bomb("unable to read input file");

  if (computeGGE(&fieldsOnBoundary, normalOutputFile, skewOutputFile, maxDerivatives, maxMultipoles, fundamental, varyDerivatives))
    return 1;

  bggexpData[0].haveData = bggexpData[1].haveData = 0;
  readBGGExpData(&bggexpData[0], normalOutputFile, "CnmS", 0);
  if (skewOutputFile)
    readBGGExpData(&bggexpData[1], skewOutputFile, "CnmC", 1);

  if (autoTuneFlags & AUTOTUNE_ACTIVE) {
    if (autoTuneFlags & AUTOTUNE_EVALONLY)
      derivatives = maxDerivatives;
    else
      derivatives = 1;
  } else
    /* no auto-tuning */
    derivatives = maxDerivatives;

  minMultipoles = 1;
  for (; derivatives <= maxDerivatives; derivatives++) {
    multipoles = minMultipoles;
    if (autoTuneFlags & AUTOTUNE_ACTIVE) {
      if (autoTuneFlags & AUTOTUNE_EVALONLY)
        multipoles = maxMultipoles;
    } else
      /* no auto-tuning */
      multipoles = maxMultipoles;

    for (; multipoles <= maxMultipoles; multipoles++) {
      if (autoTuneFlags & AUTOTUNE_ACTIVE) {
        double residual;
        if ((residual = evaluateGGEFit(&fieldsOnBoundary, &bggexpData[0],
                                       derivatives, multipoles,
                                       autoTuneSignificance, autoTuneFlags, varyDerivatives, &allResiduals)) < bestResidual) {
          bestResidual = residual;
          bestMultipoles = multipoles;
          bestDerivatives = derivatives;
          if (autoTuneFlags & AUTOTUNE_INCRONLY)
            minMultipoles = bestMultipoles;
          if (autoTuneFlags & AUTOTUNE_VERBOSE) {
            printf("New best residual of %le for m=%ld, d=%ld\n", residual, multipoles, derivatives);
            fflush(stdout);
          }
        } else {
          if (autoTuneFlags & AUTOTUNE_VERBOSE) {
            printf("Goodness of fit (%le) for m=%ld, d=%ld is not better than %le for m=%ld, d=%ld\n",
                   residual, multipoles, derivatives, bestResidual, bestMultipoles, bestDerivatives);
          }
        }
        if (autoTuneFlags & AUTOTUNE_LOG) {
          if (!SDDS_SetRowValues(&SDDS_autoTuneLog, SDDS_SET_BY_NAME | SDDS_PASS_BY_VALUE,
                                 iAutoTuneLog++,
                                 "m", multipoles * (fundamental > 0 ? fundamental : 1), "d", derivatives, "RmsError", allResiduals.rms,
                                 "MaximumError", allResiduals.max, "MadError", allResiduals.mad,
                                 "FractionalRmsError", allResiduals.fracRms,
                                 "FractionalMadError", allResiduals.fracMad,
                                 "FractionalMaximumError", allResiduals.fracMax,
                                 NULL)) {
            SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
            return (1);
          }
        }
      }
    }
  }

  freeBGGExpData(&bggexpData[0]);
  freeBGGExpData(&bggexpData[1]);

  if (autoTuneFlags & AUTOTUNE_ACTIVE && !(autoTuneFlags & AUTOTUNE_EVALONLY) &&
      computeGGE(&fieldsOnBoundary, normalOutputFile, skewOutputFile, bestDerivatives, bestMultipoles, fundamental, varyDerivatives))
    return 1;

  if (autoTuneFlags & AUTOTUNE_LOG) {
    char buffer[1024];
    snprintf(buffer, 1024, "m$bopt$n: %ld  d$bopt$n: %ld  r$bopt$n: %lg T",
             bestMultipoles, bestDerivatives, bestResidual);
    if (SDDS_SetParameters(&SDDS_autoTuneLog, SDDS_SET_BY_NAME | SDDS_PASS_BY_VALUE,
                           "OptimalMultipoles", bestMultipoles,
                           "OptimalDerivatives", bestDerivatives,
                           "OptimalResidual", bestResidual,
                           "OptimumLabel", buffer,
                           NULL) != 1) {
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
    }
    if (SDDS_WritePage(&SDDS_autoTuneLog) != 1 || SDDS_Terminate(&SDDS_autoTuneLog) != 1) {
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
      return (1);
    }
  }

  if (evaluationOutput)
    evaluateGGEAndOutput(evaluationOutput, evaluation_nRho, evaluation_nPhi,
                         normalOutputFile, skewOutputFile, &fieldsOnBoundary);

  return (0);
}

int ReadInputFile(
                  FIELDS_ON_BOUNDARY *fieldsOnBoundary,
                  char *inputFile,
                  char *zName,
                  char *phiName,
                  char *BrhoName,
                  char *BzName,
                  char *rhoName) {
  SDDS_DATASET SDDSin;
  double *z, *phi;
  int32_t rows, iphi, iz;
  double deltaPhiMin, deltaPhiMax, deltaPhi;
  double deltaZMin, deltaZMax, deltaZ;

  /* Read in Brho */
  if (SDDS_InitializeInput(&SDDSin, inputFile) != 1) {
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
    return (1);
  }
  if ((SDDS_CheckColumn(&SDDSin, BrhoName, "T", SDDS_ANY_NUMERIC_TYPE, stderr) != SDDS_CHECK_OKAY) ||
      (SDDS_CheckColumn(&SDDSin, phiName, NULL, SDDS_ANY_NUMERIC_TYPE, stderr) != SDDS_CHECK_OKAY) ||
      (SDDS_CheckColumn(&SDDSin, zName, "m", SDDS_ANY_NUMERIC_TYPE, stderr) != SDDS_CHECK_OKAY) ||
      (SDDS_CheckParameter(&SDDSin, rhoName, "m", SDDS_ANY_NUMERIC_TYPE, stderr) != SDDS_CHECK_OKAY)) {
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
    return (1);
  }
  if (BzName && strlen(BzName) &&
      SDDS_CheckColumn(&SDDSin, BzName, "T", SDDS_ANY_NUMERIC_TYPE, stderr) != SDDS_CHECK_OKAY) {
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
    return (1);
  }
  if (SDDS_ReadPage(&SDDSin) != 1) {
    fprintf(stderr, "Unable to read SDDS page\n");
    return (1);
  }
  rows = SDDS_RowCount(&SDDSin);
  if (!(fieldsOnBoundary->Brho = SDDS_GetColumnInDoubles(&SDDSin, BrhoName)) ||
      !(phi = SDDS_GetColumnInDoubles(&SDDSin, phiName)) ||
      !(z = SDDS_GetColumnInDoubles(&SDDSin, zName))) {
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
    return (1);
  }
  fieldsOnBoundary->Bz = NULL;
  if (BzName &&
      !(fieldsOnBoundary->Bz = SDDS_GetColumnInDoubles(&SDDSin, BzName))) {
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
    return (1);
  }
  if (!SDDS_GetParameterAsDouble(&SDDSin, rhoName, &(fieldsOnBoundary->rho))) {
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
    return (1);
  }

  for (iz = 1; iz < rows; iz++) {
    if (z[iz - 1] > z[iz])
      SDDS_Bomb("Data is not correctly ordered in input file. Use sddssort -column=z -column=phi on the file.");
  }

  for (iz = 1; iz < rows; iz++) {
    if (z[iz - 1] != z[iz])
      break;
  }
  fieldsOnBoundary->Nphi = iz;
  fieldsOnBoundary->Nz = rows / iz;

  if ((fieldsOnBoundary->Nz * fieldsOnBoundary->Nphi) != rows)
    SDDS_Bomb("number of z and phi values doesn't multiply to number of total rows in input file");

  for (iphi = 1; iphi < fieldsOnBoundary->Nphi; iphi++) {
    if (phi[iphi - 1] >= phi[iphi]) {
      fprintf(stderr, "Data is not correctly ordered in input file. Use sddssort -column=%s -column=%s on the file.\n",
              zName, phiName);
      return 1;
    }
  }
  for (iphi = 0; iphi < fieldsOnBoundary->Nphi; iphi++) {
    for (iz = 0; iz < fieldsOnBoundary->Nz; iz++) {
      if (phi[iz * fieldsOnBoundary->Nphi + iphi] != phi[iphi])
        SDDS_Bomb("phi data is not repeated at each z plane");
    }
  }
  for (iz = 0; iz < fieldsOnBoundary->Nz; iz++) {
    for (iphi = 0; iphi < fieldsOnBoundary->Nphi; iphi++) {
      if (z[iz * fieldsOnBoundary->Nphi + iphi] != z[iz * fieldsOnBoundary->Nphi])
        SDDS_Bomb("z data is not fixed for expected number of slots");
    }
  }

  deltaPhiMin = DBL_MAX;
  deltaPhiMax = -DBL_MAX;
  for (iphi = 1; iphi < fieldsOnBoundary->Nphi; iphi++) {
    deltaPhi = phi[iphi] - phi[iphi - 1];
    if (deltaPhi > deltaPhiMax)
      deltaPhiMax = deltaPhi;
    if (deltaPhi < deltaPhiMin)
      deltaPhiMin = deltaPhi;
  }
  if ((deltaPhiMax - deltaPhiMin) > 1e-6)
    SDDS_Bomb("delta phi varies by more than 1e-6");
  fieldsOnBoundary->dphi = (deltaPhiMin + deltaPhiMax) / 2;

  deltaZMin = DBL_MAX;
  deltaZMax = -DBL_MAX;
  for (iz = 1; iz < fieldsOnBoundary->Nz; iz++) {
    deltaZ = z[iz * fieldsOnBoundary->Nphi] - z[(iz - 1) * fieldsOnBoundary->Nphi];
    if (deltaZ > deltaZMax)
      deltaZMax = deltaZ;
    if (deltaZ < deltaZMin)
      deltaZMin = deltaZ;
  }
  if ((deltaZMax - deltaZMin) > 1e-6)
    SDDS_Bomb("delta z varies by more than 1e-6");
  fieldsOnBoundary->dz = (deltaZMin + deltaZMax) / 2;
  fieldsOnBoundary->zMin = z[0];

#ifdef DEBUG
  fprintf(stderr, "Nz = %ld, Nphi = %ld\n", (long)fieldsOnBoundary->Nz, (long)fieldsOnBoundary->Nphi);
  fprintf(stderr, "dz = %le, dphi = %le\n", fieldsOnBoundary->dz, fieldsOnBoundary->dphi);
  fprintf(stderr, "zMin = %le\n", fieldsOnBoundary->zMin);
  fprintf(stderr, "rho = %le\n", fieldsOnBoundary->rho);
#endif

  free(phi);
  free(z);

  return 0;
}

double Imp(double kR, long m) {
  double t1;

  t1 = (BesIn(fabs(kR), m - 1) + BesIn(fabs(kR), m + 1)) / 2;
  if (kR < 0)
    return t1 * ipow(-1, m + 1);
  return t1;
}

int computeGGE(
               FIELDS_ON_BOUNDARY *fieldsOnBoundary,
               char *normalOutput,
               char *skewOutput,
               long derivatives,
               long multipoles,
               long fundamental,
               long varyDerivatives) {
  COMPLEX **Brho, *Bz, **Bm, *c1;
  double *k;
  double rho, dk;
  long ik, n, iphi, iz, Nz, Nphi, offsetN, offsetS, m, iharm;
  SDDS_DATASET SDDSnormal, SDDSskew;
#ifdef DEBUG
  FILE *fp;
#endif

  Nz = fieldsOnBoundary->Nz;
  Nphi = fieldsOnBoundary->Nphi;
  rho = fieldsOnBoundary->rho;

#ifdef DEBUG
  if (fieldsOnBoundary->Bz) {
    fp = fopen("gge0.sdds", "w");
    fprintf(fp, "SDDS1\n");
    fprintf(fp, "&parameter name = phi, type = double, units = m &end\n");
    fprintf(fp, "&column name = z, type = double  &end\n");
    fprintf(fp, "&column name = Bz, type = double  &end\n");
    fprintf(fp, "&data mode = ascii, &end\n");
    for (iphi = 0; iphi < Nphi; iphi++) {
      fprintf(fp, "%le\n", PIx2 * iphi / Nphi);
      fprintf(fp, "%ld\n", Nz);
      for (iz = 0; iz < Nz; iz++)
        fprintf(fp, "%le %le\n", fieldsOnBoundary->zMin + iz * fieldsOnBoundary->dz, fieldsOnBoundary->Bz[iphi + iz * Nphi]);
    }
    fclose(fp);
  }
#endif

  /* Take FFT of Brho values vs phi for each z plane */
  Brho = malloc(sizeof(*Brho) * Nz);
#pragma omp parallel
  {
    int iz, iphi, myid;
#if !defined(NOTHREADS)
    myid = omp_get_thread_num();
#else
    myid = 0;
#endif
    for (iz = 0; iz < Nz; iz++) {
      if (iz % threads == myid) {
        Brho[iz] = malloc(sizeof(**Brho) * Nphi);
        for (iphi = 0; iphi < Nphi; iphi++) {
          Brho[iz][iphi].re = fieldsOnBoundary->Brho[iphi + iz * Nphi];
          Brho[iz][iphi].im = 0;
        }
        FFT(Brho[iz], -1, Nphi);
      }
    }
#pragma omp barrier
  }

  if (fieldsOnBoundary->Bz) {
    /* Optionally average Bz over phi for each z plane */
    Bz = malloc(sizeof(*Bz) * Nz);
    for (iz = 0; iz < Nz; iz++) {
      Bz[iz].re = Bz[iz].im = 0;
      for (iphi = 0; iphi < Nphi; iphi++)
        Bz[iz].re += fieldsOnBoundary->Bz[iphi + iz * Nphi];
      Bz[iz].re /= Nphi;
    }
  } else
    Bz = NULL;

  /* Reorganize data to give array of FFT for each frequency as a function of z, then take FFT vs z */
  Bm = malloc(sizeof(*Bm) * Nphi);
#pragma omp parallel
  {
    int iphi, iz, myid;
#if !defined(NOTHREADS)
    myid = omp_get_thread_num();
#else
    myid = 0;
#endif
    for (iphi = 0; iphi < Nphi; iphi++) {
      if (iphi % threads == myid) {
        Bm[iphi] = tmalloc(sizeof(**Bm) * Nz);
        for (iz = 0; iz < Nz; iz++) {
          Bm[iphi][iz].re = Brho[iz][iphi].re;
          Bm[iphi][iz].im = Brho[iz][iphi].im;
        }
        FFT(Bm[iphi], -1, Nz);
      }
    }
#pragma omp barrier
  }

  if (Bz)
    FFT(Bz, -1, Nz);

  /* clean up memory */
  for (iz = 0; iz < Nz; iz++)
    free(Brho[iz]);
  free(Brho);

  /* Compute k values. Upper half of the array has negative k values */
  dk = TWOPI / (fieldsOnBoundary->dz * Nz);
  k = calloc(Nz, sizeof(double));
  for (ik = 0; ik < Nz / 2; ik++)
    k[ik] = dk * ik;
  k[0] = 1e-12 * dk;
  for (ik = Nz / 2; ik < Nz; ik++)
    k[ik] = -dk * (Nz - ik);

  if (normalOutput && SetUpOutputFile(&SDDSnormal, normalOutput, 0, derivatives))
    return 1;
  if (skewOutput && SetUpOutputFile(&SDDSskew, skewOutput, 1, derivatives))
    return 1;

  c1 = malloc(Nz * sizeof(*c1));
  offsetN = 1;
  offsetS = 1;
  for (m = (Bz ? -1 : 0); m < multipoles; m++) {
    if (fundamental)
      iharm = fundamental * (2 * m + 1);
    else
      iharm = m + 1;
    if (normalOutput && m >= 0) {
      offsetN = SDDS_GetColumnIndex(&SDDSnormal, "CnmS0");
      StartPageLocal(&SDDSnormal, fieldsOnBoundary, iharm, "Normal");
    }
    if (skewOutput) {
      offsetS = SDDS_GetColumnIndex(&SDDSskew, "CnmC0");
      StartPageLocal(&SDDSskew, fieldsOnBoundary, iharm, "Skew");
    }
#ifdef DEBUG
    printf("GGE m=%ld, iharm=%ld, offsetN = %ld, offsetS = %ld\n",
           m, iharm, offsetN, offsetS);
#endif

    if (m >= 0) {
      long nLimit;
      nLimit = 2 * derivatives;
      if (varyDerivatives)
        nLimit = 2 * (derivatives - m / 2);
      for (n = 0; n < 2 * derivatives; n += 2) {
        for (ik = 0; ik < Nz; ik++) {
          double factor;
          factor = ipow(-1, n / 2) * ipow(k[ik], iharm + n - 1) / Imp(k[ik] * rho, iharm) / (ipow(2, iharm) * dfactorial(iharm)) / (0.5 * Nz * Nphi);
          c1[ik].re = -Bm[iharm][ik].im * factor;
          c1[ik].im = Bm[iharm][ik].re * factor;
        }
        FFT(c1, 1, Nz);
        for (iz = 0; iz < Nz; iz++) {
          if (normalOutput &&
              !SDDS_SetRowValues(&SDDSnormal, SDDS_SET_BY_INDEX | SDDS_PASS_BY_VALUE, iz,
                                 offsetN, n > nLimit ? 0.0 : c1[iz].re,
                                 -1)) {
            SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
            return (1);
          }
          if (skewOutput &&
              !SDDS_SetRowValues(&SDDSskew, SDDS_SET_BY_INDEX | SDDS_PASS_BY_VALUE, iz,
                                 offsetS, n > nLimit ? 0.0 : c1[iz].im,
                                 -1)) {
            SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
            return (1);
          }
        }
        offsetN++;
        offsetS++;
      }

      if (normalOutput)
        offsetN = SDDS_GetColumnIndex(&SDDSnormal, "dCnmS0/dz");
      if (skewOutput)
        offsetS = SDDS_GetColumnIndex(&SDDSskew, "dCnmC0/dz");

#ifdef DEBUG
      printf("dGGE m=%ld, offsetN = %ld, offsetS = %ld\n",
             m, offsetN, offsetS);
#endif
      for (n = 0; n < 2 * derivatives; n += 2) {
        for (ik = 0; ik < Nz; ik++) {
          double factor;
          factor = -ipow(-1, n / 2) * ipow(k[ik], iharm + n - 1) / Imp(k[ik] * rho, iharm) / (ipow(2, iharm) * dfactorial(iharm)) * k[ik] / Nz / (Nphi / 2.0);
          c1[ik].re = Bm[iharm][ik].re * factor;
          c1[ik].im = Bm[iharm][ik].im * factor;
        }
        FFT(c1, 1, Nz);
        for (iz = 0; iz < Nz; iz++) {
          if (normalOutput &&
              !SDDS_SetRowValues(&SDDSnormal, SDDS_SET_BY_INDEX | SDDS_PASS_BY_VALUE, iz,
                                 offsetN, n > nLimit ? 0.0 : c1[iz].re,
                                 -1)) {
            SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
            return (1);
          }
          if (skewOutput &&
              !SDDS_SetRowValues(&SDDSskew, SDDS_SET_BY_INDEX | SDDS_PASS_BY_VALUE, iz,
                                 offsetS, n > nLimit ? 0.0 : c1[iz].im,
                                 -1)) {
            SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
            return (1);
          }
        }
        offsetN++;
        offsetS++;
      }
    }

    if (m < 0 && skewOutput && Bz) {
      /* solenoidal terms will be included here */
      offsetS = SDDS_GetColumnIndex(&SDDSskew, "CnmC0");
      for (n = 0; n < 2 * derivatives; n += 2) {
        if ((iharm + n) > 0) {
          for (ik = 0; ik < Nz; ik++) {
            double factor;
            factor = -ipow(-1, n / 2) * ipow(k[ik], iharm + n - 1) / BesIn(k[ik] * rho, iharm) / (ipow(2, iharm) * dfactorial(iharm)) / (Nz);
            c1[ik].re = -Bz[ik].im * factor;
            c1[ik].im = Bz[ik].re * factor;
          }
          FFT(c1, 1, Nz);
          for (iz = 0; iz < Nz; iz++) {
            if (!SDDS_SetRowValues(&SDDSskew, SDDS_SET_BY_INDEX | SDDS_PASS_BY_VALUE, iz,
                                   offsetS, c1[iz].re,
                                   -1)) {
              SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
              return (1);
            }
          }
        }
        offsetS++;
      }
      offsetS = SDDS_GetColumnIndex(&SDDSskew, "dCnmC0/dz");
      for (n = 0; n < 2 * derivatives; n += 2) {
        for (ik = 0; ik < Nz; ik++) {
          double factor;
          factor = ipow(-1, n / 2) * ipow(k[ik], iharm + n) / BesIn(k[ik] * rho, iharm) / (ipow(2, iharm) * dfactorial(iharm)) / (Nz);
          c1[ik].re = Bz[ik].re * factor;
          c1[ik].im = Bz[ik].im * factor;
        }
        FFT(c1, 1, Nz);
        for (iz = 0; iz < Nz; iz++) {
          if (skewOutput &&
              !SDDS_SetRowValues(&SDDSskew, SDDS_SET_BY_INDEX | SDDS_PASS_BY_VALUE, iz,
                                 offsetS, c1[iz].re,
                                 -1)) {
            SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
            return (1);
          }
        }
        offsetS++;
      }
    }
    if (normalOutput && m >= 0 && SDDS_WritePage(&SDDSnormal) != 1) {
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
      return (1);
    }
    if (skewOutput && SDDS_WritePage(&SDDSskew) != 1) {
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
      return (1);
    }
  }

  if (normalOutput && SDDS_Terminate(&SDDSnormal) != 1) {
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
    return (1);
  }
  if (skewOutput && SDDS_Terminate(&SDDSskew) != 1) {
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
    return (1);
  }

#ifdef DEBUG
  printf("Finished writing GGE output files\n");
  fflush(stdout);
#endif

  if (Bz)
    free(Bz);
  for (iphi = 0; iphi < Nphi; iphi++)
    free(Bm[iphi]);
  free(Bm);

  return (0);
}

void readFieldMap(char *fieldMapFile, FIELD_MAP *fmData) {
  SDDS_DATASET SDDSin;

  if (!SDDS_InitializeInputFromSearchPath(&SDDSin, fieldMapFile))
    SDDS_Bomb("unable to read field input file");
  if (SDDS_CheckColumn(&SDDSin, "Bx", "T", SDDS_ANY_NUMERIC_TYPE, stderr) != SDDS_CHECK_OKAY ||
      SDDS_CheckColumn(&SDDSin, "By", "T", SDDS_ANY_NUMERIC_TYPE, stderr) != SDDS_CHECK_OKAY ||
      SDDS_CheckColumn(&SDDSin, "Bz", "T", SDDS_ANY_NUMERIC_TYPE, stderr) != SDDS_CHECK_OKAY)
    SDDS_Bomb("Didn't find required field columns Bx, By, Bz in T");
  if (SDDS_CheckColumn(&SDDSin, "x", "m", SDDS_ANY_NUMERIC_TYPE, stderr) != SDDS_CHECK_OKAY ||
      SDDS_CheckColumn(&SDDSin, "y", "m", SDDS_ANY_NUMERIC_TYPE, stderr) != SDDS_CHECK_OKAY ||
      SDDS_CheckColumn(&SDDSin, "z", "m", SDDS_ANY_NUMERIC_TYPE, stderr) != SDDS_CHECK_OKAY)
    SDDS_Bomb("Didn't find required coordinate columns x, y, z in T");
  if (SDDS_ReadPage(&SDDSin) <= 0 ||
      !(fmData->x = SDDS_GetColumnInDoubles(&SDDSin, "x")) ||
      !(fmData->y = SDDS_GetColumnInDoubles(&SDDSin, "y")) ||
      !(fmData->z = SDDS_GetColumnInDoubles(&SDDSin, "z")) ||
      !(fmData->Bx = SDDS_GetColumnInDoubles(&SDDSin, "Bx")) ||
      !(fmData->By = SDDS_GetColumnInDoubles(&SDDSin, "By")) ||
      !(fmData->Bz = SDDS_GetColumnInDoubles(&SDDSin, "Bz")))
    SDDS_Bomb("unable to get data from field input file");
  if (!(fmData->n = SDDS_CountRowsOfInterest(&SDDSin)) || fmData->n < 1)
    SDDS_Bomb("field map file has insufficient data");
  SDDS_Terminate(&SDDSin);
}

void freeBGGExpData(BGGEXP_DATA *bggexpData) {
  long im, id;
  if (bggexpData->haveData) {
    for (im = 0; im < bggexpData->nm; im++) {
      for (id = 0; id < bggexpData->nGradients; id++) {
        if (bggexpData->Cmn[im][id])
          free(bggexpData->Cmn[im][id]);
        if (bggexpData->dCmn_dz[im][id])
          free(bggexpData->dCmn_dz[im][id]);
        bggexpData->Cmn[im][id] = bggexpData->dCmn_dz[im][id] = NULL;
      }
    }
    free(bggexpData->m);
    bggexpData->m = NULL;
    bggexpData->nm = bggexpData->nGradients = 0;
  }
}

#define BUFSIZE 1024

void readBGGExpData(BGGEXP_DATA *bggexpData, char *filename, char *nameFragment, short skew) {
  SDDS_DATASET SDDSin;
  char buffer[BUFSIZE];
  long im, ic, nc, readCode, nz;
  int32_t m;
  short xCenterPresent = 0, yCenterPresent = 0, xMaxPresent = 0, yMaxPresent = 0;

  bggexpData->haveData = 1;

  if (!SDDS_InitializeInput(&SDDSin, filename)) {
    fprintf(stderr, "Unable to read file %s\n", filename);
    exit(1);
  }

  /* Check presence of z column */
  if (SDDS_CheckColumn(&SDDSin, "z", "m", SDDS_ANY_FLOATING_TYPE, stderr) != SDDS_CHECK_OK) {
    fprintf(stderr, "Unable to find floating-point column \"z\" with units \"m\" in file %s\n", filename);
    exit(1);
  }

  /* Check presence of Cnm* columns */
  ic = 0;
  while (1) {
    snprintf(buffer, BUFSIZE, "%s%ld", nameFragment, 2 * ic);
    if (SDDS_CheckColumn(&SDDSin, buffer, NULL, SDDS_ANY_FLOATING_TYPE, NULL) != SDDS_CHECK_OK)
      break;
    ic++;
  }
  if (ic == 0) {
    fprintf(stderr, "Unable to find any floating-point columns %s* in file %s\n",
            nameFragment, filename);
    exit(1);
  }
  nc = ic;

  /* Check for presence of matching dCnmXXX/dz columns */
  for (ic = 0; ic < nc; ic++) {
    snprintf(buffer, BUFSIZE, "d%s%ld/dz", nameFragment, 2 * ic);
    if (SDDS_CheckColumn(&SDDSin, buffer, NULL, SDDS_ANY_FLOATING_TYPE, stderr) != SDDS_CHECK_OK)
      break;
  }
  if (ic != nc) {
    fprintf(stderr, "Unable to find matching floating-point columns dCnm*/dz in file %s\n", filename);
    exit(1);
  }

  bggexpData->nm = bggexpData->nz = 0;
  bggexpData->nGradients = nc;
  bggexpData->m = NULL;
  bggexpData->Cmn = NULL;
  bggexpData->dCmn_dz = NULL;

  im = nz = 0;
  bggexpData->zMin = DBL_MAX;
  bggexpData->zMax = -DBL_MAX;
  bggexpData->xCenter = bggexpData->yCenter = 0;
  bggexpData->xMax = bggexpData->yMax = -1;
  while ((readCode = SDDS_ReadPage(&SDDSin)) > 0) {
    if (!SDDS_GetParameterAsLong(&SDDSin, "m", &m) || (m < 1 && !skew) || (m < 0 && skew)) {
      fprintf(stderr, "Problem with value of m (m<%d) for page %ld of file %s\n",
              (skew ? 0 : 1), readCode, filename);
      exit(1);
    }
    if (readCode == 1) {
      long iz;
      double dz0, dz, *z, zMin, zMax;
      if ((xCenterPresent && !SDDS_GetParameterAsDouble(&SDDSin, "xCenter", &(bggexpData->xCenter))) ||
          (yCenterPresent && !SDDS_GetParameterAsDouble(&SDDSin, "yCenter", &(bggexpData->yCenter)))) {
        fprintf(stderr, "Problem getting xCenter or yCenter values from file %s\n",
                filename);
        exit(1);
      }
      if ((xMaxPresent && !SDDS_GetParameterAsDouble(&SDDSin, "xMax", &(bggexpData->xMax))) ||
          (yMaxPresent && !SDDS_GetParameterAsDouble(&SDDSin, "yMax", &(bggexpData->yMax)))) {
        fprintf(stderr, "Problem getting xMax or yMax values from file %s\n",
                filename);
        exit(1);
      }
      if ((nz = SDDS_RowCount(&SDDSin)) <= 1) {
        fprintf(stderr, "Too few z values in file %s\n", filename);
        exit(1);
      }
      if (!(z = SDDS_GetColumnInDoubles(&SDDSin, "z"))) {
        fprintf(stderr, "Problem reading column z from %s\n", filename);
        exit(1);
      }
      find_min_max(&zMin, &zMax, z, nz);
      if (zMin < bggexpData->zMin)
        bggexpData->zMin = zMin;
      if (zMax > bggexpData->zMax)
        bggexpData->zMax = zMax;
      dz0 = z[1] - z[0];
      for (iz = 1; iz < nz; iz++) {
        dz = z[iz] - z[iz - 1];
        if (dz <= 0 || fabs(dz0 / dz - 1) > 1e-6) {
          fprintf(stderr, "Data not uniformly and monotonically increasing in z column from %s\n", filename);
          exit(1);
        }
      }
      free(z);
      bggexpData->dz = dz0;
    } else {
      if (nz != SDDS_RowCount(&SDDSin)) {
        fprintf(stderr, "Inconsistent number of z values in file %s\n", filename);
        exit(1);
      }
    }
    if (!(bggexpData->m = SDDS_Realloc(bggexpData->m,
                                       sizeof(*bggexpData->m) * (im + 1))) ||
        !(bggexpData->Cmn = SDDS_Realloc(bggexpData->Cmn,
                                         sizeof(*bggexpData->Cmn) * (im + 1))) ||
        !(bggexpData->dCmn_dz = SDDS_Realloc(bggexpData->dCmn_dz,
                                             sizeof(*bggexpData->dCmn_dz) * (im + 1)))) {
      fprintf(stderr, "Memory allocation failure (1) loading data from file %s\n", filename);
      exit(1);
    }
    bggexpData->Cmn[im] = NULL;
    bggexpData->dCmn_dz[im] = NULL;
    if (!(bggexpData->Cmn[im] = malloc(sizeof(*bggexpData->Cmn[im]) * nc)) ||
        !(bggexpData->dCmn_dz[im] = malloc(sizeof(*bggexpData->dCmn_dz[im]) * nc)))
      fprintf(stderr, "Memory allocation failure (2) loading data from file %s\n", filename);

    for (ic = 0; ic < nc; ic++) {
      snprintf(buffer, BUFSIZE, "%s%ld", nameFragment, 2 * ic);
      if (!(bggexpData->Cmn[im][ic] = SDDS_GetColumnInDoubles(&SDDSin, buffer))) {
        SDDS_PrintErrors(stdout, SDDS_VERBOSE_PrintErrors);
        fprintf(stderr, "Problem reading column %s from %s\n", buffer, filename);
        exit(1);
      }
      snprintf(buffer, BUFSIZE, "d%s%ld/dz", nameFragment, 2 * ic);
      if (!(bggexpData->dCmn_dz[im][ic] = SDDS_GetColumnInDoubles(&SDDSin, buffer))) {
        SDDS_PrintErrors(stdout, SDDS_VERBOSE_PrintErrors);
        fprintf(stderr, "Problem reading column %s from %s\n", buffer, filename);
        exit(1);
      }
    }

    bggexpData->m[im] = m;
    im++;
  }

  SDDS_Terminate(&SDDSin);

  bggexpData->nm = im;
  bggexpData->nz = nz;
}

double evaluateGGEFit(
                      FIELDS_ON_BOUNDARY *fieldsOnBoundary,
                      BGGEXP_DATA bggexpData[2],
                      long derivatives,
                      long multipoles,
                      double significance,
                      unsigned long flags,
                      long varyDerivatives,
                      ALL_RESIDUALS *allResiduals) {
  double *residualWorst, *residualSum, *residualSum2, *maxField;
  double r;
  int i;

  r = fieldsOnBoundary->rho;
  residualWorst = calloc(threads, sizeof(*residualWorst));
  residualSum = calloc(threads, sizeof(*residualSum));
  residualSum2 = calloc(threads, sizeof(*residualSum2));
  maxField = calloc(threads, sizeof(*maxField));

#pragma omp parallel
  {
    long ns, ig, m, im;
    long iphi, iz, myid;
    double Br, Bphi, Bz, field;
    double phi, residual;

#if !defined(NOTHREADS)
    myid = omp_get_thread_num();
#else
    myid = 0;
#endif

    residualWorst[myid] = residualSum[myid] = residualSum2[myid] = 0;
    maxField[myid] = -1;

    for (iz = 0; iz < fieldsOnBoundary->Nz; iz++) {
      if (iz % threads != myid)
        continue;
      for (iphi = 0; iphi < fieldsOnBoundary->Nphi; iphi++) {
        phi = (TWOPI * iphi) / fieldsOnBoundary->Nphi;

        /* Compute fields */
        Br = Bphi = Bz = 0;
        for (ns = 0; ns < 2; ns++) {
          /* ns=0 => normal, ns=1 => skew */
          if (!bggexpData[ns].haveData)
            continue;

          for (im = 0; im < bggexpData[ns].nm && im < multipoles; im++) {
            double mfact, term, sin_mphi, cos_mphi;
            long ndLimit;
            m = bggexpData[ns].m[im];
            ndLimit = MIN(bggexpData[ns].nGradients, derivatives);
            if (varyDerivatives)
              ndLimit -= m / 2;
            mfact = dfactorial(m);
            sin_mphi = sin(m * phi);
            cos_mphi = cos(m * phi);
            if (ns == 0) {
              /* normal */
              for (ig = 0; ig < ndLimit; ig++) {
                term = ipow(-1, ig) * mfact / (ipow(2, 2 * ig) * factorial(ig) * factorial(ig + m)) * ipow(r, 2 * ig + m - 1);
                Bz += term * bggexpData[ns].dCmn_dz[im][ig][iz] * r * sin_mphi;
                term *= bggexpData[ns].Cmn[im][ig][iz];
                Br += term * (2 * ig + m) * sin_mphi;
                Bphi += m * term * cos_mphi;
              }
            } else {
              /* skew */
              if (m == 0) {
                Bz += bggexpData[ns].dCmn_dz[im][0][iz]; // on-axis Bz from m=ig=0 term
                for (ig = 1; ig < ndLimit; ig++) {
                  term = ipow(-1, ig) * mfact / (ipow(2, 2 * ig) * factorial(ig) * factorial(ig + m)) * ipow(r, 2 * ig + m - 1);
                  Bz += term * bggexpData[ns].dCmn_dz[im][ig][iz] * r;
                  Br += term * (2 * ig + m) * bggexpData[ns].Cmn[im][ig][iz];
                }
              } else {
                for (ig = 0; ig < ndLimit; ig++) {
                  term = ipow(-1, ig) * mfact / (ipow(2, 2 * ig) * factorial(ig) * factorial(ig + m)) * ipow(r, 2 * ig + m - 1);
                  Bz += term * bggexpData[ns].dCmn_dz[im][ig][iz] * r * cos_mphi;
                  term *= bggexpData[ns].Cmn[im][ig][iz];
                  Br += term * (2 * ig + m) * cos_mphi;
                  Bphi -= m * term * sin_mphi;
                }
              }
            }
          }
        }
        field = fabs(Br);
        if (field > maxField[myid])
          maxField[myid] = field;
        residual = fabs(Br - fieldsOnBoundary->Brho[iz * fieldsOnBoundary->Nphi + iphi]);
        if (residual > residualWorst[myid])
          residualWorst[myid] = residual;
        residualSum2[myid] += sqr(residual);
        residualSum[myid] += residual;
      }
    }
#pragma omp barrier
  }

  for (i = 1; i < threads; i++) {
    residualWorst[0] = MAX(residualWorst[0], residualWorst[i]);
    residualSum[0] += residualSum[i];
    residualSum2[0] += residualSum2[i];
    maxField[0] = MAX(maxField[0], maxField[i]);
  }

  allResiduals->rms = sqrt(residualSum2[0] / (fieldsOnBoundary->Nz * fieldsOnBoundary->Nphi));
  allResiduals->mad = residualSum[0] / (fieldsOnBoundary->Nz * fieldsOnBoundary->Nphi);
  allResiduals->max = residualWorst[0];
  if (maxField[0] > 0) {
    allResiduals->fracRms = allResiduals->rms / maxField[0];
    allResiduals->fracMad = allResiduals->mad / maxField[0];
    allResiduals->fracMax = allResiduals->max / maxField[0];
  } else
    allResiduals->fracRms =
      allResiduals->fracMad =
      allResiduals->fracMax = -DBL_MAX;

  if (flags & AUTOTUNE_RMS) {
    residualWorst[0] = allResiduals->rms;
  } else if (flags & AUTOTUNE_MAV) {
    residualWorst[0] = allResiduals->mad;
  }

  return residualWorst[0] > significance ? residualWorst[0] : 0.0;
}

int evaluateGGEAndOutput(char *outputFile, long nrho, long nphi,
                         char *normalFile, char *skewFile, FIELDS_ON_BOUNDARY *fieldsOnBoundary) {
  double B[3], Br, Bphi;
  double z, r, phi;
  long ns, ig, m, im;
  BGGEXP_DATA bggexpData[2];
  SDDS_DATASET SDDSout;
  long iphi, iz, irow, irho;

  bggexpData[0].haveData = bggexpData[1].haveData = 0;
  readBGGExpData(&bggexpData[0], normalFile, "CnmS", 0);
  if (skewFile)
    readBGGExpData(&bggexpData[1], skewFile, "CnmC", 1);

  if (nrho < 0)
    nrho = 1;
  if (nphi < 0)
    nphi = fieldsOnBoundary->Nphi;

  if (SDDS_InitializeOutput(&SDDSout, SDDS_BINARY, 1, NULL, NULL, outputFile) != 1 ||
      !SDDS_DefineSimpleColumn(&SDDSout, "x", "m", SDDS_DOUBLE) ||
      !SDDS_DefineSimpleColumn(&SDDSout, "y", "m", SDDS_DOUBLE) ||
      !SDDS_DefineSimpleColumn(&SDDSout, "z", "m", SDDS_DOUBLE) ||
      !SDDS_DefineSimpleColumn(&SDDSout, "rho", "m", SDDS_DOUBLE) ||
      !SDDS_DefineSimpleColumn(&SDDSout, "phi", NULL, SDDS_DOUBLE) ||
      !SDDS_DefineSimpleColumn(&SDDSout, "Bx", "T", SDDS_DOUBLE) ||
      !SDDS_DefineSimpleColumn(&SDDSout, "By", "T", SDDS_DOUBLE) ||
      !SDDS_DefineSimpleColumn(&SDDSout, "Bz", "T", SDDS_DOUBLE) ||
      !SDDS_DefineSimpleColumn(&SDDSout, "Bphi", "T", SDDS_DOUBLE) ||
      !SDDS_DefineSimpleColumn(&SDDSout, "Brho", "T", SDDS_DOUBLE) ||
      !SDDS_WriteLayout(&SDDSout) ||
      !SDDS_StartPage(&SDDSout, fieldsOnBoundary->Nz * nrho * nphi)) {
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
    return 1;
  }

  irow = 0;
  for (iz = 0; iz < fieldsOnBoundary->Nz; iz++) {
    z = fieldsOnBoundary->zMin + fieldsOnBoundary->dz * iz;
    for (irho = 0; irho < nrho; irho++) {
      r = (irho + 1) * (fieldsOnBoundary->rho / nrho);
      for (iphi = 0; iphi < nphi; iphi++) {
        phi = (TWOPI * iphi) / nphi;

        /* Compute fields */
        Br = Bphi = B[0] = B[1] = B[2] = 0;
        for (ns = 0; ns < 2; ns++) {
          /* ns=0 => normal, ns=1 => skew */
          if (!bggexpData[ns].haveData)
            continue;

          for (im = 0; im < bggexpData[ns].nm; im++) {
            double mfact, term, sin_mphi, cos_mphi;
            m = bggexpData[ns].m[im];
            mfact = dfactorial(m);
            sin_mphi = sin(m * phi);
            cos_mphi = cos(m * phi);
            if (ns == 0) {
              /* normal */
              for (ig = 0; ig < bggexpData[ns].nGradients; ig++) {
                term = ipow(-1, ig) * mfact / (ipow(2, 2 * ig) * factorial(ig) * factorial(ig + m)) * ipow(r, 2 * ig + m - 1);
                B[2] += term * bggexpData[ns].dCmn_dz[im][ig][iz] * r * sin_mphi;
                term *= bggexpData[ns].Cmn[im][ig][iz];
                Br += term * (2 * ig + m) * sin_mphi;
                Bphi += m * term * cos_mphi;
              }
            } else {
              /* skew */
              if (m == 0) {
                B[2] += bggexpData[ns].dCmn_dz[im][0][iz]; // on-axis Bz from m=ig=0 term
                for (ig = 1; ig < bggexpData[ns].nGradients; ig++) {
                  term = ipow(-1, ig) * mfact / (ipow(2, 2 * ig) * factorial(ig) * factorial(ig + m)) * ipow(r, 2 * ig + m - 1);
                  B[2] += term * bggexpData[ns].dCmn_dz[im][ig][iz] * r;
                  Br += term * (2 * ig + m) * bggexpData[ns].Cmn[im][ig][iz];
                }
              } else {
                for (ig = 0; ig < bggexpData[ns].nGradients; ig++) {
                  term = ipow(-1, ig) * mfact / (ipow(2, 2 * ig) * factorial(ig) * factorial(ig + m)) * ipow(r, 2 * ig + m - 1);
                  B[2] += term * bggexpData[ns].dCmn_dz[im][ig][iz] * r * cos_mphi;
                  term *= bggexpData[ns].Cmn[im][ig][iz];
                  Br += term * (2 * ig + m) * cos_mphi;
                  Bphi -= m * term * sin_mphi;
                }
              }
            }
          }
        }
        B[0] = Br * cos(phi) - Bphi * sin(phi);
        B[1] = Br * sin(phi) + Bphi * cos(phi);
        if (SDDS_SetRowValues(&SDDSout, SDDS_SET_BY_NAME | SDDS_PASS_BY_VALUE, irow++,
                              "x", r * cos(phi), "y", r * sin(phi), "z", z,
                              "phi", phi, "rho", r,
                              "Bx", B[0], "By", B[1], "Bz", B[2],
                              "Bphi", Bphi, "Brho", Br,
                              NULL) != 1) {
          SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
          return (1);
        }
      }
    }
  }
  if (!SDDS_WritePage(&SDDSout) || !SDDS_Terminate(&SDDSout)) {
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
    return (1);
  }

  freeBGGExpData(&bggexpData[0]);
  freeBGGExpData(&bggexpData[1]);
  return 0;
}

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
      if (((long)labs(order)) % 2)
        result = -gsl_sf_bessel_In(order, -x);
      else
        result = gsl_sf_bessel_In(order, -x);
    }
  }
  return result;
}

int SetUpOutputFile(SDDS_DATASET *SDDSout, char *filename, long skew, long derivatives) {
  long n;
  char tag[2] = "SC";
  char name[1024], units[1024];

  if (skew)
    skew = 1;

  if (SDDS_InitializeOutput(SDDSout, SDDS_BINARY, 1, NULL,
                            skew ? "computeCBGGE skew output" : "computeCBGGE normal output", filename) != 1) {
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
    return (1);
  }
  if ((SDDS_DefineSimpleParameter(SDDSout, "m", NULL, SDDS_LONG) != 1) ||
      (SDDS_DefineSimpleParameter(SDDSout, "rho", "m", SDDS_DOUBLE) != 1) ||
      (SDDS_DefineSimpleColumn(SDDSout, "z", "m", SDDS_DOUBLE) != 1)) {
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
    return (1);
  }

  for (n = 0; n < 2 * derivatives; n += 2) {
    sprintf(name, "Cnm%c%ld", tag[skew], n);
    if ((2 * n - 1) < 0)
      sprintf(units, "T/m$a(m-%ld)$n", -(2 * n - 1));
    else if ((2 * n - 1) == 0)
      sprintf(units, "T/m$am$n");
    else
      sprintf(units, "T/m$a(m+%ld)$n", (2 * n - 1));
    if (SDDS_DefineSimpleColumn(SDDSout, name, units, SDDS_DOUBLE) != 1) {
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
      return (1);
    }
  }
  for (n = 0; n < 2 * derivatives; n += 2) {
    sprintf(name, "dCnm%c%ld/dz", tag[skew], n);
    if ((2 * n - 2) < 0)
      sprintf(units, "T/m$a(m-%ld)$n", -(2 * n - 2));
    else if ((2 * n - 2) == 0)
      sprintf(units, "T/m$am$n");
    else
      sprintf(units, "T/m$a(m+%ld)$n", (2 * n - 2));
    if (SDDS_DefineSimpleColumn(SDDSout, name, units, SDDS_DOUBLE) != 1) {
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
      return (1);
    }
  }
  if (!SDDS_WriteLayout(SDDSout)) {
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
    return (1);
  }
  return 0;
}

int StartPageLocal(SDDS_DATASET *SDDSout,
              FIELDS_ON_BOUNDARY *fob,
              long m,
              char *type) {
  long iz, index;
  if (!SDDS_StartPage(SDDSout, fob->Nz)) {
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
    return (1);
  }
  if (SDDS_SetParameters(SDDSout, SDDS_SET_BY_NAME | SDDS_PASS_BY_VALUE,
                         "m", m, "rho", fob->rho,
                         NULL) != 1) {
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
    return (1);
  }
  index = SDDS_GetColumnIndex(SDDSout, "z");
#ifdef DEBUG
  printf("StartPageLocal %s: dz = %le, Nz = %ld, column index = %ld\n",
         type, fob->dz, fob->Nz, index);
#endif
  for (iz = 0; iz < fob->Nz; iz++) {
    if (!SDDS_SetRowValues(SDDSout, SDDS_SET_BY_INDEX | SDDS_PASS_BY_VALUE,
                           iz, index, iz * fob->dz + fob->zMin, -1)) {
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
      return (1);
    }
  }
  return 0;
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
