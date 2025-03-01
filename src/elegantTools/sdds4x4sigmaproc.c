/*************************************************************************\
 * Copyright (c) 2018 The University of Chicago, as Operator of Argonne
 * National Laboratory.
 * Copyright (c) 2018 The Regents of the University of California, as
 * Operator of Los Alamos National Laboratory.
 * This file is distributed subject to a Software License Agreement found
 * in the file LICENSE that is included with this distribution.
\*************************************************************************/

/* program: sdds4x4sigmaproc
   purpose: process elegant runs to determine the 4x4 sigma matrix and thus emittances of beam,
   assuming that the transport line being varied has no coupling.
   Also processes experimental data using elegant-computed matrices.
   Michael Borland, 2018.

   Uses the following equations:
   S11(e) = R11^2*S11(s) + 2*R11*R12*S12(s) + R12^2*S22(s)
   S33(e) = R33^2*S33(s) + 2*R33*R34*S34(s) + R34^2*S44(s)
   S13(e) = R11*R33*S13(s) + R11*R34*S14(s) + R12*R33*S23(s) + R12*R34*S24(s)
   where (e) is the exit (observation  point) and (s) is the start.
   We assume that the transport matrix Rij is block-diagonal (i.e., no cross-plane terms).
   Thus, all the information about coupling terms in the input beam comes from the
   measurement of the x-y correlation at the observation point.

   There are N measurements of Sij(e) as the matrix Rij is varied.
   To solve for Sij(s), we organize the data into the form
   M = P.F
   The equation is solved by finding F = Inverse(Transpose(P) P) Transpose(P) M

   Take the example of the x plane, S11(e) = R11^2*S11(s) + 2*R11*R12*S12(s) + R12^2*S22(s), first.

   M is a N x 1 column vector with the measurements.
   Tranpose(M) = ( S11_1(e) ... S11_N(e) )

   P is a N x 3 matrix that contains the R matrix data from simulation
   P = ( R11_1^2  2*R11_1*R12_1 R12_1^2
   .
   .
   .
   R11_N^2  2*R11_N*R12_N R12_N^2 )

   F is a 3 x 1 matrix of the initial sigma matrix elements
   Transpose(F) = (S11(s) S12(s) S22(s))

   For the y plane, just replace R11 by R33, S11 by S33, etc.

   For the x-y terms,
   Transpose(M) = (S13_1(e) ... S13_N(e))
   P is N x 4
   P = ( R11_1*R33_1 R11_1*R34_1 R12_1*R33_1 R12_1*R34_1
   .
   .
   .
   R11_N*R33_N R11_N*R34_N R12_N*R33_N R12_N*R34_N)
   F is 4 x 1
   Transpose(F) = (S13(s) S14(s) S23(s) S24(s))

*/

#include "mdb.h"
#include "matlib.h"
#include "match_string.h"
#include "SDDS.h"
#include "scan.h"
#include "table.h"
#include <time.h>
#include <memory.h>

#define SET_N_ERROR_SETS 0
#define SET_SEED 1
#define SET_VERBOSITY 2
#define SET_DEVIATION_LIMIT 3
#define SET_RESOLUTION 4
#define SET_LIMIT_MODE 5
#define SET_PIPE 6
#define N_OPTIONS 7

char *option[N_OPTIONS] = {
  "nerrorsets",
  "seed",
  "verbosity",
  "deviationlimit",
  "resolution",
  "limitmode",
  "pipe",
};

#define USAGE "sddsemitproc\n\
 [<inputfile>] [<outputfile>] [-pipe=[input][,output]]\n\
 [-nErrorSets=<number> [-seed=<integer>]]\n\
 [-deviationLimit=<valueInSigma>] [-limitMode={resolution | zero}[{,reject}]\n\
 [-resolution=<xResolutionm>,<yResolutionm>]\n\
 [-seed=integer] [-verbosity=level]\n\n\
Program by Michael Borland. (This is version 1, August 2018)"

static char *additional_help[] = {
  USAGE,
  "\nThis program computes the 4x4 sigma matrix from simulated or experimental",
  "beamsize data, using simulated data for transport matrices for a beamline that",
  "itself has no coupling terms.",
  "inputfile is an elegant \"final\" parameters output file (SDDS format),",
  "which contains the R matrix and simulated beam sigmas as a function",
  "of some variable or variables.  It may also contain experimental data",
  "that has been obtained separately and put into the elegant output file",
  "(e.g., using sddsxref). \n",
  "-nErrorsSets is used to specify the number of randomizations of the data",
  "    to use in estimating errors in the sigma-matrix determination.",
  "    Error levels are given by columns S11Sigma, S33Sigma, and S13Sigma, representing the",
  "    error in determination of the S11, S33, and S13 values for a given row.",
  "-limitMode is used to specify what to do if a randomized measurement lies",
  "    too far from the fit.",
  "-deviationLimit is used to define what \"too far\" from the fit means.",
  "-resolution allows specification of the measurement resolution,",
  "    which is subtracted in quadrature from the sigma.",
  NULL};

#define N_INCREMENT 10

double solve_normal_form(MATRIX *F, MATRIX *sF, MATRIX *P, MATRIX *M, MATRIX *C, double *s2_fit);
double solve_normal_form_opt(MATRIX *F, MATRIX *sF, MATRIX *P, MATRIX *M, MATRIX *C, double dev_limit,
                             long *n_used, double *s2_fit);
double propagate_errors_for_emittance(double **Sigma, double **Covar);
void set_up_covariance_matrix(MATRIX *K, double *sigma, double *uncert, long n_configs, long equal_weights);
double estimate_uncertainty(double *uncert, MATRIX *S, MATRIX *sS, MATRIX *R, MATRIX *s2, MATRIX *K,
                            double dev_limit, long n_configs, double uncert_min, double *fit_sig2_return);
long SetSigmaData(SDDS_DATASET *SDDSout, char *dataName, MATRIX *s2, char *fitName, double *fitSigSqr,
                  long configs);
double *solveForSigmaMatrix(int i, int j, double *SMeasured, double *SSigma,
                            double *R11, double *R12, double *R33, double *R34,
                            long nConfigs, long nErrorSets, long *nValidReturn, double resolution,
                            double SijSum[4][4], double SijSum2[4][4]);

#define GAUSSIAN_ERRORS 0
#define UNIFORM_ERRORS 1
#define N_ERROR_TYPES 2
char *error_type[N_ERROR_TYPES] = {
  "gaussian", "uniform"};

#define LIMIT_AT_RESOLUTION 0
#define LIMIT_AT_ZERO 1
#define N_LIMIT_OPTIONS 2
char *limit_option[N_LIMIT_OPTIONS] = {
  "resolution", "zero"};

#define MAX_N_TRIES 100

int main(
         int argc,
         char **argv) {
  SDDS_TABLE SDDSin, SDDSout;
  double *R11, *R12 = NULL, *R33 = NULL, *R34 = NULL;
  long n_configs;
  double *S11, *S11Sigma;               /* measured <x^2> for ith config, plus measurement uncertainty */
  double *S33 = NULL, *S33Sigma = NULL; /* measured <y^2> for ith config, plus measurement uncertainty */
  double *S13 = NULL, *S13Sigma = NULL; /* measured <xy> for ith config, plus measurement uncertainty */
  double SijSum[4][4], SijSum2[4][4];   /* sum of inferred sigma matrix values, their squares */
  double *S11Fit, *S33Fit, *S13Fit;
  long seed;
  long nErrorSets, nGoodFits[3];
  SCANNED_ARG *scanned;
  long i_arg, i;
  char *input, *output;
  double x_resol, y_resol;
  double deviationLimit = 0;
  // double x_limit=0.0, y_limit=0.0;
  long limit_code;
  // long reject_at_limit;
  long verbosity;
  unsigned long pipeFlags;

  argc = scanargs(&scanned, argc, argv);
  if (argc < 2 || argc > (2 + N_OPTIONS)) {
    for (i = 0;; i++) {
      if (!additional_help[i])
        break;
      puts(additional_help[i]);
    }
    exit(1);
  }

  input = output = NULL;
  nErrorSets = 0;
  seed = -1;
  x_resol = y_resol = 0;
  limit_code = LIMIT_AT_ZERO;
  // reject_at_limit = 0;
  verbosity = 0;
  pipeFlags = 0;

  for (i_arg = 1; i_arg < argc; i_arg++) {
    if (scanned[i_arg].arg_type == OPTION) {
      switch (match_string(scanned[i_arg].list[0], option,
                           N_OPTIONS, 0)) {
        /* process options here */
      case SET_N_ERROR_SETS:
        if (scanned[i_arg].n_items != 2 ||
            !sscanf(scanned[i_arg].list[1], "%ld", &nErrorSets) ||
            nErrorSets < 0)
          bomb("invalid -nErrorSets syntax", USAGE);
        break;
      case SET_SEED:
        if (scanned[i_arg].n_items != 2 ||
            !sscanf(scanned[i_arg].list[1], "%ld", &seed) ||
            seed < 0)
          bomb("invalid -seed syntax", USAGE);
        break;
      case SET_DEVIATION_LIMIT:
        if (scanned[i_arg].n_items != 2 ||
            !sscanf(scanned[i_arg].list[1], "%lf", &deviationLimit) ||
            deviationLimit <= 0)
          bomb("invalid -deviationLimit syntax", USAGE);
        break;
      case SET_RESOLUTION:
        if (scanned[i_arg].n_items != 3 ||
            !sscanf(scanned[i_arg].list[1], "%lf", &x_resol) ||
            !sscanf(scanned[i_arg].list[2], "%lf", &y_resol) ||
            x_resol < 0 || y_resol < 0)
          bomb("invalid -resolution syntax", USAGE);
        break;
      case SET_LIMIT_MODE:
        if (scanned[i_arg].n_items < 2 || scanned[i_arg].n_items > 4 ||
            (limit_code = match_string(scanned[i_arg].list[1], limit_option,
                                       N_LIMIT_OPTIONS, 0)) < 0)
          bomb("invalid -limit_mode syntax", USAGE);
        if (scanned[i_arg].n_items == 3) {
          if (scanned[i_arg].list[2][0] == 'r') {
            // reject_at_limit = 1;
          } else
            bomb("invalid -limit_mode syntax", USAGE);
        }
        break;
      case SET_VERBOSITY:
        if (scanned[i_arg].n_items != 2 ||
            !sscanf(scanned[i_arg].list[1], "%ld", &verbosity) ||
            verbosity < 0)
          bomb("invalid -verbosity syntax", USAGE);
        break;
      case SET_PIPE:
        if (!processPipeOption(scanned[i_arg].list + 1, scanned[i_arg].n_items - 1, &pipeFlags))
          SDDS_Bomb("invalid -pipe syntax");
        break;
      default:
        bomb("unknown option given", USAGE);
        break;
      }
    } else {
      if (!input)
        input = scanned[i_arg].list[0];
      else if (!output)
        output = scanned[i_arg].list[0];
      else
        bomb("too many filenames given", USAGE);
    }
  }

  processFilenames("sdds4x4sigmaproc", &input, &output, pipeFlags, 0, NULL);

  // x_limit = (limit_code==LIMIT_AT_ZERO?0.0:x_resol);
  // y_limit = (limit_code==LIMIT_AT_ZERO?0.0:y_resol);

  if (seed < 0)
    /* generate seed from system clock */
    seed = (int)time(NULL);
  random_1(-seed);

  if (nErrorSets <= 2)
    SDDS_Bomb("number of error sets must be >2");

  if (!SDDS_InitializeInput(&SDDSin, input))
    SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors | SDDS_VERBOSE_PrintErrors);
  if (!SDDS_InitializeOutput(&SDDSout, SDDS_BINARY, 1, NULL, NULL, output))
    SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors | SDDS_VERBOSE_PrintErrors);
  if (!SDDS_DefineSimpleParameter(&SDDSout, "S11", "m$a2$n", SDDS_DOUBLE) ||
      !SDDS_DefineSimpleParameter(&SDDSout, "S12", "m", SDDS_DOUBLE) ||
      !SDDS_DefineSimpleParameter(&SDDSout, "S22", "", SDDS_DOUBLE) ||
      !SDDS_DefineSimpleParameter(&SDDSout, "S33", "m$a2$n", SDDS_DOUBLE) ||
      !SDDS_DefineSimpleParameter(&SDDSout, "S34", "m", SDDS_DOUBLE) ||
      !SDDS_DefineSimpleParameter(&SDDSout, "S44", "", SDDS_DOUBLE) ||
      !SDDS_DefineSimpleParameter(&SDDSout, "S13", "m$a2$n", SDDS_DOUBLE) ||
      !SDDS_DefineSimpleParameter(&SDDSout, "S14", "m", SDDS_DOUBLE) ||
      !SDDS_DefineSimpleParameter(&SDDSout, "S23", "m", SDDS_DOUBLE) ||
      !SDDS_DefineSimpleParameter(&SDDSout, "S24", "", SDDS_DOUBLE) ||
      !SDDS_DefineSimpleParameter(&SDDSout, "goodFits", "", SDDS_LONG) ||
      !SDDS_DefineSimpleParameter(&SDDSout, "averageFitPoints", "", SDDS_DOUBLE) ||
      !SDDS_DefineSimpleColumn(&SDDSout, "S11Data", "m$a2$n", SDDS_DOUBLE) ||
      !SDDS_DefineSimpleColumn(&SDDSout, "S11Fit", "m$a2$n", SDDS_DOUBLE) ||
      !SDDS_DefineSimpleColumn(&SDDSout, "S33Data", "m$a2$n", SDDS_DOUBLE) ||
      !SDDS_DefineSimpleColumn(&SDDSout, "S33Fit", "m$a2$n", SDDS_DOUBLE) ||
      !SDDS_DefineSimpleColumn(&SDDSout, "S13Data", "m$a2$n", SDDS_DOUBLE) ||
      !SDDS_DefineSimpleColumn(&SDDSout, "S13Fit", "m$a2$n", SDDS_DOUBLE)) {
    SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors | SDDS_VERBOSE_PrintErrors);
  }

  if (SDDS_GetColumnIndex(&SDDSin, "S11") < 0 ||
      SDDS_GetColumnIndex(&SDDSin, "S33") < 0 ||
      SDDS_GetColumnIndex(&SDDSin, "S13") < 0 ||
      SDDS_GetColumnIndex(&SDDSin, "R11") < 0 ||
      SDDS_GetColumnIndex(&SDDSin, "R12") < 0 ||
      SDDS_GetColumnIndex(&SDDSin, "R33") < 0 ||
      SDDS_GetColumnIndex(&SDDSin, "R34") < 0)
    SDDS_Bomb("input file missing required quantities. Need S11, S33, S13, R11, R12, R33, and R34.");

  if (SDDS_GetColumnIndex(&SDDSin, "S11Sigma") < 0 ||
      SDDS_GetColumnIndex(&SDDSin, "S33Sigma") < 0 ||
      SDDS_GetColumnIndex(&SDDSin, "S13Sigma") < 0)
    SDDS_Bomb("input file missing required quantities. Need S11Sigma, S33Sigma, and S13Sigma.");

  if (!SDDS_WriteLayout(&SDDSout))
    SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors | SDDS_VERBOSE_PrintErrors);

  while (SDDS_ReadTable(&SDDSin) > 0) {
    n_configs = SDDS_CountRowsOfInterest(&SDDSin);
    if (!(R11 = SDDS_GetColumn(&SDDSin, "R11")) ||
        !(R12 = SDDS_GetColumn(&SDDSin, "R12")) ||
        !(R33 = SDDS_GetColumn(&SDDSin, "R33")) ||
        !(R34 = SDDS_GetColumn(&SDDSin, "R34")))
      SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors | SDDS_VERBOSE_PrintErrors);

    if (!(S11 = SDDS_GetColumn(&SDDSin, "S11")) ||
        !(S33 = SDDS_GetColumn(&SDDSin, "S33")) ||
        !(S13 = SDDS_GetColumn(&SDDSin, "S13")))
      SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors | SDDS_VERBOSE_PrintErrors);

    if (!(S11Sigma = SDDS_GetColumn(&SDDSin, "S11Sigma")) ||
        !(S33Sigma = SDDS_GetColumn(&SDDSin, "S33Sigma")) ||
        !(S13Sigma = SDDS_GetColumn(&SDDSin, "S13Sigma")))
      SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors | SDDS_VERBOSE_PrintErrors);

    if (n_configs < 4)
      continue;

    if (!SDDS_StartPage(&SDDSout, n_configs) || !SDDS_CopyColumns(&SDDSout, &SDDSin))
      SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors | SDDS_VERBOSE_PrintErrors);

    S11Fit = solveForSigmaMatrix(1, 1, S11, S11Sigma, R11, R12, R33, R34, n_configs,
                                 nErrorSets, &nGoodFits[0], x_resol, SijSum, SijSum2);
    S33Fit = solveForSigmaMatrix(3, 3, S33, S33Sigma, R11, R12, R33, R34, n_configs,
                                 nErrorSets, &nGoodFits[1], y_resol, SijSum, SijSum2);
    S13Fit = solveForSigmaMatrix(1, 3, S13, S13Sigma, R11, R12, R33, R34, n_configs,
                                 nErrorSets, &nGoodFits[2], sqrt(x_resol * y_resol), SijSum, SijSum2);

    if (nGoodFits[0] && nGoodFits[1] && nGoodFits[2]) {
      if (!SDDS_StartPage(&SDDSout, n_configs))
        SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors | SDDS_VERBOSE_PrintErrors);

      if (!SDDS_SetParameters(&SDDSout, SDDS_SET_BY_NAME | SDDS_PASS_BY_VALUE,
                              "S11", SijSum[0][0] / nGoodFits[0], "S12", SijSum[0][1] / nGoodFits[0], "S22", SijSum[1][1] / nGoodFits[0],
                              NULL))
        SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors | SDDS_VERBOSE_PrintErrors);

      if (!SDDS_SetParameters(&SDDSout, SDDS_SET_BY_NAME | SDDS_PASS_BY_VALUE,
                              "S33", SijSum[2][2] / nGoodFits[1], "S34", SijSum[2][3] / nGoodFits[1], "S44", SijSum[3][3] / nGoodFits[1],
                              NULL))
        SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors | SDDS_VERBOSE_PrintErrors);

      if (!SDDS_SetParameters(&SDDSout, SDDS_SET_BY_NAME | SDDS_PASS_BY_VALUE,
                              "S13", SijSum[0][2] / nGoodFits[2], "S14", SijSum[0][3] / nGoodFits[2],
                              "S23", SijSum[1][2] / nGoodFits[2], "S24", SijSum[1][3] / nGoodFits[2],
                              NULL))
        SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors | SDDS_VERBOSE_PrintErrors);

      if (!SDDS_SetColumn(&SDDSout, SDDS_SET_BY_NAME, S11, n_configs, "S11Data") ||
          !SDDS_SetColumn(&SDDSout, SDDS_SET_BY_NAME, S33, n_configs, "S33Data") ||
          !SDDS_SetColumn(&SDDSout, SDDS_SET_BY_NAME, S13, n_configs, "S13Data") ||
          !SDDS_SetColumn(&SDDSout, SDDS_SET_BY_NAME, S11Fit, n_configs, "S11Fit") ||
          !SDDS_SetColumn(&SDDSout, SDDS_SET_BY_NAME, S33Fit, n_configs, "S33Fit") ||
          !SDDS_SetColumn(&SDDSout, SDDS_SET_BY_NAME, S13Fit, n_configs, "S13Fit") ||
          !SDDS_WritePage(&SDDSout))
        SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors | SDDS_VERBOSE_PrintErrors);
    }
    free(S11Fit);
    free(S33Fit);
    free(S13Fit);
  }

  if (!SDDS_Terminate(&SDDSin) || !SDDS_Terminate(&SDDSout))
    SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors | SDDS_VERBOSE_PrintErrors);

  return (0);
}

double *solveForSigmaMatrix(int i, int j,
                            double *SMeasured, double *SSigma, /* measured Sij, sigma of same */
                            double *R11, double *R12, double *R33, double *R34,
                            long nConfigs,
                            long nErrorSets,
                            long *nValidReturn,
                            double resolution,
                            double SijSum[4][4],
                            double SijSum2[4][4]) {
  /* Sets up problem to solve M = P.F for one of the three cases described in the comment at the top
     of this file.
  */
  long iConfig, iErrorSet, nUnknowns = 0, nConfigUsed, nValid;
  MATRIX *P, *F, *M, *K, *FSigma;
  double *s2Fit, *s2FitSaved;
  long i1, i2;
  short badPoint;

  if (i == j && (i == 1 || i == 3)) {
    nUnknowns = 3;
  } else if (i == 1 && j == 3) {
    nUnknowns = 4;
  } else
    bomb("Invalid parameters for solveForSigmaMatrix. Seek expert help.", NULL);

  m_alloc(&P, nConfigs, nUnknowns);
  m_alloc(&F, nUnknowns, 1);
  m_alloc(&FSigma, nUnknowns, nUnknowns);
  m_alloc(&M, nConfigs, 1);
  m_alloc(&K, nConfigs, nConfigs);

  /* Set up P matrix */
  for (iConfig = 0; iConfig < nConfigs; iConfig++) {
    if (i == j) {
      if (i == 1) {
        P->a[iConfig][0] = sqr(R11[iConfig]);
        P->a[iConfig][1] = 2 * R11[iConfig] * R12[iConfig];
        P->a[iConfig][2] = sqr(R12[iConfig]);
        for (i1 = 0; i1 < 2; i1++)
          for (i2 = 0; i2 < 2; i2++)
            SijSum[i1][i2] = SijSum2[i1][i2] = 0;
      } else {
        P->a[iConfig][0] = sqr(R33[iConfig]);
        P->a[iConfig][1] = 2 * R33[iConfig] * R34[iConfig];
        P->a[iConfig][2] = sqr(R34[iConfig]);
        for (i1 = 2; i1 < 4; i1++)
          for (i2 = 2; i2 < 4; i2++)
            SijSum[i1][i2] = SijSum2[i1][i2] = 0;
      }
    } else {
      P->a[iConfig][0] = R11[iConfig] * R33[iConfig];
      P->a[iConfig][1] = R11[iConfig] * R34[iConfig];
      P->a[iConfig][2] = R12[iConfig] * R33[iConfig];
      P->a[iConfig][3] = R12[iConfig] * R34[iConfig];
      for (i1 = 0; i1 < 2; i1++)
        for (i2 = 2; i2 < 4; i2++)
          SijSum[i1][i2] = SijSum2[i1][i2] =
            SijSum[i2][i1] = SijSum2[i2][i1] = 0;
    }
  }

  /* Setup inverse covariance matrix */
  m_zero(K);
  for (iConfig = 0; iConfig < nConfigs; iConfig++)
    K->a[iConfig][iConfig] = 1 / SSigma[iConfig];

  s2Fit = tmalloc(sizeof(*s2Fit) * nConfigs);
  s2FitSaved = NULL;
  nValid = 0;
  for (iErrorSet = 0; iErrorSet < nErrorSets; iErrorSet++) {
    /* Set up problem for new error instance */
    for (iConfig = badPoint = 0; iConfig < nConfigs; iConfig++)
      if ((M->a[iConfig][0] = SMeasured[iConfig] + gauss_rn_lim(0.0, SSigma[iConfig], 2, random_1) - sqr(resolution)) <= 0 && i == j)
        badPoint = 1;
    if (badPoint)
      continue;
    solve_normal_form_opt(F, FSigma, P, M, K, 0.0, &nConfigUsed, s2Fit);
    if (nConfigUsed > 4) {
      nValid++;
      s2FitSaved = s2Fit;
      if (i == j) {
        if (i == 1) {
          SijSum[0][0] += F->a[0][0];
          SijSum[0][1] += F->a[1][0];
          SijSum[1][1] += F->a[2][0];
          SijSum2[0][0] += sqr(F->a[0][0]);
          SijSum2[0][1] += sqr(F->a[1][0]);
          SijSum2[1][1] += sqr(F->a[2][0]);
        } else {
          SijSum[2][2] += F->a[0][0];
          SijSum[2][3] += F->a[1][0];
          SijSum[3][3] += F->a[2][0];
          SijSum2[2][2] += sqr(F->a[0][0]);
          SijSum2[2][3] += sqr(F->a[1][0]);
          SijSum2[3][3] += sqr(F->a[2][0]);
        }
      } else {
        SijSum[0][2] += F->a[0][0];
        SijSum[0][3] += F->a[1][0];
        SijSum[1][2] += F->a[2][0];
        SijSum[1][3] += F->a[3][0];
        SijSum2[0][2] += sqr(F->a[0][0]);
        SijSum2[0][3] += sqr(F->a[1][0]);
        SijSum2[1][2] += sqr(F->a[2][0]);
        SijSum2[1][3] += sqr(F->a[3][0]);
      }
    }
  }

  m_free(&P);
  m_free(&F);
  m_free(&FSigma);
  m_free(&M);
  m_free(&K);

  *nValidReturn = nValid;
  return s2FitSaved;
}

double solve_normal_form_opt(
                             MATRIX *F,  /* Mx1 matrix of Fit coefficients (returned) */
                             MATRIX *sF, /* MxM matrix of errors in Fit coefficients (returned) */
                             MATRIX *P,  /* NxM matrix of Parameters of fit (provided) */
                             MATRIX *M,  /* Nx1 column matrix of Measured quantities (provided) */
                             /* M = P.F is the equation being solved */
                             MATRIX *K,        /* NxN inverse covariance matrix for Measured quantities.
                                                  K[i][j] = delta(i,j)/uncert[i]^2 */
                             double dev_limit, /* limit on deviation for any point used in final fit */
                             long *n_used,     /* number of points used in final fit (returned) */
                             double *s2_fit    /* sigma^2 from fit (returned) */
                             ) {
  long i, j, i_good;
  double rms_error, error;
  double *s2_fit2;
  long *index, n_good, *good;
  MATRIX *Pc, *Mc, *Kc;

  rms_error = solve_normal_form(F, sF, P, M, K, s2_fit);

  if (dev_limit == 0) {
    *n_used = M->n;
    return (rms_error);
  }

  good = tmalloc(sizeof(*good) * M->n);
  if (dev_limit < 0)
    dev_limit = -dev_limit * rms_error;
  for (i = n_good = 0; i < M->n; i++) {
    error = fabs(sqrt(s2_fit[i]) - sqrt(M->a[i][0]));
    if (dev_limit > error) {
      n_good++;
      good[i] = 1;
    } else {
      s2_fit[i] = -1; /* used to mark excluded points for caller */
      good[i] = 0;
    }
  }
  if (n_good == 0) {
    *n_used = 0;
    free(good);
    return (0.0);
  }

  m_alloc(&Pc, n_good, P->m);
  m_alloc(&Mc, n_good, 1);
  m_alloc(&Kc, n_good, n_good);
  m_zero(Kc);
  index = tmalloc(sizeof(int) * n_good);
  s2_fit2 = tmalloc(sizeof(*s2_fit2) * n_good);

  for (i = i_good = 0; i < M->n; i++) {
    if (good[i]) {
      index[i_good] = i;
      for (j = 0; j < P->m; j++)
        Pc->a[i_good][j] = P->a[i][j];
      for (j = 0; j < Mc->m; j++)
        Mc->a[i_good][j] = M->a[i][j];
      Kc->a[i_good][i_good] = K->a[i][i];
      i_good++;
    }
  }

  *n_used = n_good;
  rms_error = solve_normal_form(F, sF, Pc, Mc, Kc, s2_fit2);
  for (i = 0; i < n_good; i++)
    s2_fit[index[i]] = s2_fit2[i];

  free(good);
  free(s2_fit2);
  free(index);
  m_free(&Pc);
  m_free(&Mc);
  m_free(&Kc);

  return (rms_error);
}

double solve_normal_form(
                         MATRIX *F,  /* Mx1 matrix of Fit coefficients (returned) */
                         MATRIX *sF, /* MxM covariance matrix for Fit coefficients (returned) */
                         MATRIX *P,  /* NxM matrix of Parameters of fit (provided) */
                         MATRIX *M,  /* Nx1 column matrix of Measured sigma^2 (provided) */
                         MATRIX *K,  /* NxN inverse covariance matrix for Measured quantities (provided) .
                                        K[i][j] = delta(i,j)/uncert[i]^2 */
                         /* M = P.F is the equation being solved                          */
                         /* the solution is F = inverse(transpose(P) K P) tranpose(P) K M = T.M */
                         /*           and  sF = T.K.transpose(T) */
                         double *s2_fit /* sigma^2 from fit (returned) */
                         ) {
  long n, m, i;
  MATRIX *Pt, *Pt_K, *Pt_K_P, *Inv_Pt_K_P, *Inv_PtKP_PtK;
  MATRIX *Mp, *Tt, *TC, *C;
  double error, rms_error;

  n = M->n; /* n is the number of experimental measurements */
  m = P->m; /* m is 3 or 6, the number of unknowns to determine */
  m_alloc(&Pt, m, n);
  m_alloc(&Pt_K, m, n);
  m_alloc(&Pt_K_P, m, m);
  m_alloc(&Inv_Pt_K_P, m, m);
  m_alloc(&Inv_PtKP_PtK, m, n);
  m_alloc(&Mp, M->n, M->m);
  m_alloc(&Tt, n, m);
  m_alloc(&TC, m, n);
  m_alloc(&C, K->n, K->m);

  /* find the fit */
  if (!m_trans(Pt, P)) {
    fprintf(stderr, "matrix error--call was: m_trans(Pt, P)\n");
    return -1;
  }
  if (!m_mult(Pt_K, Pt, K)) {
    fprintf(stderr, "matrix error--call was: m_mult(Pt_K, Pt, K)\n");
    return -1;
  }
  if (!m_mult(Pt_K_P, Pt_K, P)) {
    fprintf(stderr, "matrix error--call was: m_mult(Pt_K_P, Pt_K, P)\n");
    return -1;
  }
  if (!m_invert(Inv_Pt_K_P, Pt_K_P)) {
    fprintf(stderr, "matrix error--call was: m_invert(Inv_Pt_K_P, Pt_K_P)\n");
    return -1;
  }

  if (!m_mult(Inv_PtKP_PtK, Inv_Pt_K_P, Pt_K)) {
    fprintf(stderr, "matrix error--call was: m_mult(Inv_PtKP_PtK, Inv_Pt_K_P, Pt_K)\n");
    return -1;
  }
  if (!m_mult(F, Inv_PtKP_PtK, M)) {
    fprintf(stderr, "matrix error--call was: m_mult(F, Inv_PtKP_PtK, M)\n");
    return -1;
  }
  m_zero(sF);
  if (m_invert(C, K)) {
    if (!m_trans(Tt, Inv_PtKP_PtK)) {
      fprintf(stderr, "matrix error--call was: m_trans(Tt, Inv_PtKP_PtK)\n");
      return -1;
    }
    if (!m_mult(TC, Inv_PtKP_PtK, C)) {
      fprintf(stderr, "matrix error--call was: m_mult(TC, Inv_PtKP_PtK, C)\n");
      return -1;
    }
    if (!m_mult(sF, TC, Tt)) {
      fprintf(stderr, "matrix error--call was: m_mult(sF, TC, Tt)\n");
      fprintf(stderr, "sF: %d x %d\n", sF->n, sF->m);
      fprintf(stderr, "TC: %d x %d\n", TC->n, TC->m);
      fprintf(stderr, "Tt: %d x %d\n", Tt->n, Tt->m);
      return -1;
    }
  }

  /* evaluate the fit */
  if (!m_mult(Mp, P, F)) {
    fprintf(stderr, "matrix error--call was: m_mult(Mp, P, F)\n");
    return -1;
  }
  for (i = rms_error = 0; i < Mp->n; i++) {
    s2_fit[i] = Mp->a[i][0];
    error = fabs(Mp->a[i][0]) - fabs(M->a[i][0]);
    rms_error += sqr(error);
  }
  if (Mp->n > 3)
    rms_error = sqrt(rms_error / (Mp->n - 3));

  m_free(&Pt);
  m_free(&Pt_K);
  m_free(&Pt_K_P);
  m_free(&Inv_Pt_K_P);
  m_free(&Inv_PtKP_PtK);
  m_free(&Mp);
  m_free(&Tt);
  m_free(&TC);
  m_free(&C);

  return (rms_error);
}

void set_up_covariance_matrix(MATRIX *K, double *sigma, double *uncert, long n_configs, long equal_weights)
/* actually, K is the inverse of the covariance matrix */
{
  long i_config;

  for (i_config = 0; i_config < n_configs; i_config++) {
    if (sigma[i_config] == 0 || uncert[i_config] == 0 || equal_weights)
      K->a[i_config][i_config] = 1;
    else
      K->a[i_config][i_config] = 1. / sqr(2 * sigma[i_config] * uncert[i_config]);
  }
}

double estimate_uncertainty(double *uncert, MATRIX *S, MATRIX *sS, MATRIX *R, MATRIX *s2,
                            MATRIX *K, double dev_limit, long n_configs, double uncert_min, double *fit_sig2_return) {
  double md, *fit_sig2;
  long n_used, i_config;

  if (fit_sig2_return)
    fit_sig2 = fit_sig2_return;
  else
    fit_sig2 = tmalloc(sizeof(*fit_sig2) * n_configs);

  /* find initial fit with supplied covariance matrix */
  md = solve_normal_form_opt(S, sS, R, s2, K, 0.0, &n_used, fit_sig2);
  if (!md || !n_used)
    bomb("unable to find initial fit (1)", NULL);

  /* calculate new covariance matrix */
  for (i_config = 0; i_config < n_configs; i_config++)
    K->a[i_config][i_config] = 1. / sqr(2 * md) / s2->a[i_config][0];

  /* do second fit, excluding points that lie to far out */
  md = solve_normal_form_opt(S, sS, R, s2, K, dev_limit, &n_used, fit_sig2);
  if (!md || !n_used)
    bomb("unable to find initial fit (2)", NULL);

  /* calculate new covariance matrix */
  if (uncert_min && md < uncert_min)
    md = uncert_min;
  for (i_config = 0; i_config < n_configs; i_config++) {
    K->a[i_config][i_config] = 1. / sqr(2 * md) / s2->a[i_config][0];
    uncert[i_config] = md;
  }

  if (!fit_sig2_return)
    free(fit_sig2);

  return (md);
}

long SetSigmaData(SDDS_DATASET *SDDSout, char *dataName, MATRIX *s2,
                  char *fitName, double *fitSigSqr, long configs) {
  static double *buffer = NULL;
  long i;

  if (!(buffer = SDDS_Realloc(buffer, configs * sizeof(*buffer))))
    return 0;

  for (i = 0; i < configs; i++) {
    if (s2->a[i][0] < 0)
      buffer[i] = 0;
    else
      buffer[i] = sqrt(s2->a[i][0]);
  }
  if (!SDDS_SetColumn(SDDSout, SDDS_SET_BY_NAME, buffer, configs, dataName))
    return 0;

  for (i = 0; i < configs; i++) {
    if (fitSigSqr[i] < 0)
      buffer[i] = 0;
    else
      buffer[i] = sqrt(fitSigSqr[i]);
  }
  if (!SDDS_SetColumn(SDDSout, SDDS_SET_BY_NAME, buffer, configs, fitName))
    return 0;

  return 1;
}
