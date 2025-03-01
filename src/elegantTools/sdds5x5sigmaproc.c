/*************************************************************************\
 * Copyright (c) 2018 The University of Chicago, as Operator of Argonne
 * National Laboratory.
 * Copyright (c) 2018 The Regents of the University of California, as
 * Operator of Los Alamos National Laboratory.
 * This file is distributed subject to a Software License Agreement found
 * in the file LICENSE that is included with this distribution.
\*************************************************************************/

/* program: sdds5x5sigmaproc
   purpose: process elegant runs to determine the 5x5 sigma matrix and thus emittances of beam,
   assuming that the transport line being varied has no coupling.
   Also processes experimental data using elegant-computed matrices.
   Michael Borland, 2018.

   We assume that the measurement system does not contain vertical bending, so R36=R46=0.
   Uses the following equations:
   S11(e) = R11^2*S11(s) + 2*R11*R12*S12(s) + R12^2*S22(s) + 2*R11*R16*S16(s) + 2*R12*R16*S26(s) + R16^2*S66(s)
   S33(e) = R33^2*S33(s) + 2*R33*R34*S34(s) + R34^2*S44(s)
   S13(e) = R11*R33*S13(s) + R11*R34*S14(s) + R12*R33*S23(s) + R12*R34*S24(s) +  R16*R33*S36(s) + R16*R34*S46(s)
   where (e) is the exit (observation  point) and (s) is the start.
   We assume that the transport matrix Rij is block-diagonal (i.e., no cross-plane terms).
   Thus, all the information about coupling terms in the input beam comes from the
   measurement of the x-y correlation at the observation point.

   There are N measurements of Sij(e) as the matrix Rij is varied.
   To solve for Sij(s), we organize the data into the form
   M = P.F
   The equation is solved by finding F = Inverse(Transpose(P) P) Transpose(P) M

   Take the example of the x plane, first.

   M is a N x 1 column vector with the measurements.
   Tranpose(M) = ( S11_1(e) ... S11_N(e) )

   P is a N x 6 matrix that contains the R matrix data from simulation
   P = ( R11_1^2  2*R11_1*R12_1 R12_1^2  2*R11_1*R16_1  2*R12_1*R16_1  R16_1^2
   .
   .
   .
   R11_N^2  2*R11_N*R12_N R12_N^2  2*R11_N*R16_N  2*R12_N*R16_N  R16_N^2)

   F is a 3 x 1 matrix of the initial sigma matrix elements
   Transpose(F) = (S11(s) S12(s) S22(s) S16(s) S26(s) S66(s))

   For the y plane, just replace R11 by R33, S11 by S33, etc. and eliminate the last three columns of P and
   the last three rows of F.

   For the x-y terms,
   Transpose(M) = (S13_1(e) ... S13_N(e))
   P is N x 6
   P = ( R11_1*R33_1 R11_1*R34_1 R12_1*R33_1 R12_1*R34_1 R16_1*R33_1 R16_1*R34_1
   .
   .
   .
   R11_N*R33_N R11_N*R34_N R12_N*R33_N R12_N*R34_N R16_N*R33_N R16_N*R34_N)
   F is 6 x 1
   Transpose(F) = (S13(s) S14(s) S23(s) S24(s) S36(s) S46(s))

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
#define SET_ENERGY_SPREAD 6
#define SET_PIPE 7
#define SET_SCALE 8
#define N_OPTIONS 9

char *option[N_OPTIONS] = {
  "nerrorsets",
  "seed",
  "verbosity",
  "deviationlimit",
  "resolution",
  "limitmode",
  "energyspread",
  "pipe",
  "scale",
};

#define USAGE "sdds5x5sigmaproc\n\
 [<inputfile>] [<outputfile>] [-pipe=[input][,output]]\n\
 [-nErrorSets=<number> [-seed=<integer>]]\n\
 [-resolution=<xResolutionm>,<yResolutionm>] [-energySpread=<fractionalValue>]\n\
 [-scale=<xFactor>,<yFactor>] [-seed=integer] [-verbosity=level]\n\n\
Program by Michael Borland. (This is version 4, July 2022.)"

static char *additional_help[] = {
  USAGE,
  "\nThis program computes the 5x5 sigma matrix from simulated or experimental",
  "beamsize data, using simulated data for transport matrices for a beamline that",
  "itself has no coupling terms.",
  "inputfile is an elegant \"final\" parameters output file (SDDS format),",
  "which contains the R matrix and simulated beam sigmas as a function",
  "of some variable or variables.  It may also contain experimental data",
  "that has been obtained separately and put into the elegant output file",
  "(e.g., using sddsxref). \n",
  "-nErrorsSets is used to specify the number of randomizations of the data",
  "    to use in estimating errors in the sigma-matrix determination.",
  "    Error levels are given by columns S11StDev, S33StDev, and S13StDev, giving"
  "    the standard deviation of the measured S11, S33, and S13 values for a given row.",
  "    If set to zero (the default), then errors are ignored.",
  "-limitMode is used to specify what to do if a randomized measurement lies",
  "    too far from the fit.",
  "-deviationLimit is used to define what \"too far\" from the fit means.",
  "-resolution allows specification of the measurement resolution,",
  "    which is subtracted in quadrature from the sigma.",
  "-scale allows scaling the measured beam sizes by the given factors,",
  "    which may help identify calibration errors.",
  "-energySpread allows fixing the value of the fractional energy spread,",
  "    which may help get better results.",
  NULL};

#define N_INCREMENT 10

double solve_normal_form(MATRIX *F, MATRIX *sF, MATRIX *P, MATRIX *M, MATRIX *C, double *s2_fit);
double propagate_errors_for_emittance(double **Sigma, double **Covar);
void set_up_covariance_matrix(MATRIX *K, double *sigma, double *uncert, long nConfigs, long equal_weights);
double estimate_uncertainty(double *uncert, MATRIX *S, MATRIX *sS, MATRIX *R, MATRIX *s2, MATRIX *K,
                            double dev_limit, long nConfigs, double uncert_min, double *fit_sig2_return);
long SetSigmaData(SDDS_DATASET *SDDSout, char *dataName, MATRIX *s2, char *fitName, double *fitSigSqr,
                  long configs);
void solveForSigmaMatrix(int i, int j, double *SMeasured, double *SStDev,
                         double *R11, double *R12, double *R33, double *R34, double *R16, double *R26,
                         long nConfigs, double Sij[5][5], double *SFit, double energySpread);

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
  double *R11, *R12, *R33, *R34, *R16, *R26;
  double *S11, *S11StDev;     /* measured <x^2> for ith config, plus measurement uncertainty */
  double *S33, *S33StDev;     /* measured <y^2> for ith config, plus measurement uncertainty */
  double *S13, *S13StDev;     /* measured <xy> for ith config, plus measurement uncertainty */
  double *S11e, *S33e, *S13e; /* with errors added for Monte Carlo */
  double Sij[5][5];
  double SijSum[5][5], SijSum2[5][5]; /* sum of inferred sigma matrix values, their squares */
  double SbetaSum[4][4], emitSum[2], betaSum[2], alphaSum[2], etaSum[2], etapSum[2];
  double SbetaSum2[4][4], emitSum2[2], betaSum2[2], alphaSum2[2], etaSum2[2], etapSum2[2];
  double *S11Fit, *S33Fit, *S13Fit;
  double *S11FitSum, *S33FitSum, *S13FitSum;
  double *S11FitSum2, *S33FitSum2, *S13FitSum2;
  long seed, nErrorSets;
  SCANNED_ARG *scanned;
  long i_arg, i, j, iConfig, iSet, nConfigs;
  char *input, *output;
  double x_resol, y_resol, x_scale, y_scale;
  double deviationLimit = 0;
  long verbosity;
  unsigned long pipeFlags;
  double energySpread = 0;

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
  verbosity = 0;
  pipeFlags = 0;
  x_scale = y_scale = 1;

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
      case SET_VERBOSITY:
        if (scanned[i_arg].n_items != 2 ||
            !sscanf(scanned[i_arg].list[1], "%ld", &verbosity) ||
            verbosity < 0)
          bomb("invalid -verbosity syntax", USAGE);
        break;
      case SET_ENERGY_SPREAD:
        if (scanned[i_arg].n_items != 2 ||
            !sscanf(scanned[i_arg].list[1], "%lf", &energySpread) ||
            energySpread <= 0)
          bomb("invalid -energySpread syntax", USAGE);
        break;
      case SET_SCALE:
        if (scanned[i_arg].n_items != 3 ||
            !sscanf(scanned[i_arg].list[1], "%lf", &x_scale) ||
            !sscanf(scanned[i_arg].list[2], "%lf", &y_scale) ||
            x_scale <= 0 || y_scale <= 0)
          bomb("invalid -scale syntax", USAGE);
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

  processFilenames("sdds5x5sigmaproc", &input, &output, pipeFlags, 0, NULL);

  if (seed < 0)
    /* generate seed from system clock */
    seed = (int)time(NULL);
  random_1(-seed);

  if (nErrorSets > 0) {
    if (nErrorSets <= 2)
      SDDS_Bomb("number of error sets must be >2 if nonzero");
  } else
    nErrorSets = 1; /* no errors */

  if (!SDDS_InitializeInput(&SDDSin, input))
    SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors | SDDS_VERBOSE_PrintErrors);
  if (!SDDS_InitializeOutput(&SDDSout, SDDS_BINARY, 1, NULL, NULL, output))
    SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors | SDDS_VERBOSE_PrintErrors);
  if (!SDDS_DefineSimpleParameter(&SDDSout, "S11", "m$a2$n", SDDS_DOUBLE) ||
      !SDDS_DefineSimpleParameter(&SDDSout, "S12", "m", SDDS_DOUBLE) ||
      !SDDS_DefineSimpleParameter(&SDDSout, "S22", "", SDDS_DOUBLE) ||
      !SDDS_DefineSimpleParameter(&SDDSout, "S66", "", SDDS_DOUBLE) ||
      !SDDS_DefineSimpleParameter(&SDDSout, "S33", "m$a2$n", SDDS_DOUBLE) ||
      !SDDS_DefineSimpleParameter(&SDDSout, "S34", "m", SDDS_DOUBLE) ||
      !SDDS_DefineSimpleParameter(&SDDSout, "S44", "", SDDS_DOUBLE) ||
      !SDDS_DefineSimpleParameter(&SDDSout, "S13", "m$a2$n", SDDS_DOUBLE) ||
      !SDDS_DefineSimpleParameter(&SDDSout, "S14", "m", SDDS_DOUBLE) ||
      !SDDS_DefineSimpleParameter(&SDDSout, "S23", "m", SDDS_DOUBLE) ||
      !SDDS_DefineSimpleParameter(&SDDSout, "S24", "", SDDS_DOUBLE) ||
      !SDDS_DefineSimpleParameter(&SDDSout, "S16", "m", SDDS_DOUBLE) ||
      !SDDS_DefineSimpleParameter(&SDDSout, "S26", "", SDDS_DOUBLE) ||
      !SDDS_DefineSimpleParameter(&SDDSout, "S36", "m", SDDS_DOUBLE) ||
      !SDDS_DefineSimpleParameter(&SDDSout, "S46", "", SDDS_DOUBLE) ||
      !SDDS_DefineSimpleParameter(&SDDSout, "ex", "m", SDDS_DOUBLE) ||
      !SDDS_DefineSimpleParameter(&SDDSout, "ey", "m", SDDS_DOUBLE) ||
      !SDDS_DefineSimpleParameter(&SDDSout, "Sdelta", "", SDDS_DOUBLE) ||
      !SDDS_DefineSimpleParameter(&SDDSout, "betax", "m", SDDS_DOUBLE) ||
      !SDDS_DefineSimpleParameter(&SDDSout, "betay", "m", SDDS_DOUBLE) ||
      !SDDS_DefineSimpleParameter(&SDDSout, "alphax", "", SDDS_DOUBLE) ||
      !SDDS_DefineSimpleParameter(&SDDSout, "alphay", "", SDDS_DOUBLE) ||
      !SDDS_DefineSimpleParameter(&SDDSout, "etax", "m", SDDS_DOUBLE) ||
      !SDDS_DefineSimpleParameter(&SDDSout, "etay", "m", SDDS_DOUBLE) ||
      !SDDS_DefineSimpleParameter(&SDDSout, "etaxp", "", SDDS_DOUBLE) ||
      !SDDS_DefineSimpleParameter(&SDDSout, "etayp", "", SDDS_DOUBLE) ||
      !SDDS_DefineSimpleColumn(&SDDSout, "S11Data", "m$a2$n", SDDS_DOUBLE) ||
      !SDDS_DefineSimpleColumn(&SDDSout, "S11Fit", "m$a2$n", SDDS_DOUBLE) ||
      !SDDS_DefineSimpleColumn(&SDDSout, "S33Data", "m$a2$n", SDDS_DOUBLE) ||
      !SDDS_DefineSimpleColumn(&SDDSout, "S33Fit", "m$a2$n", SDDS_DOUBLE) ||
      !SDDS_DefineSimpleColumn(&SDDSout, "S13Data", "m$a2$n", SDDS_DOUBLE) ||
      !SDDS_DefineSimpleColumn(&SDDSout, "S13Fit", "m$a2$n", SDDS_DOUBLE)) {
    SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors | SDDS_VERBOSE_PrintErrors);
  }

  if (nErrorSets > 2 &&
      (!SDDS_DefineSimpleParameter(&SDDSout, "S11Sigma", "m$a2$n", SDDS_DOUBLE) ||
       !SDDS_DefineSimpleParameter(&SDDSout, "S12Sigma", "m", SDDS_DOUBLE) ||
       !SDDS_DefineSimpleParameter(&SDDSout, "S22Sigma", "", SDDS_DOUBLE) ||
       !SDDS_DefineSimpleParameter(&SDDSout, "S66Sigma", "", SDDS_DOUBLE) ||
       !SDDS_DefineSimpleParameter(&SDDSout, "S33Sigma", "m$a2$n", SDDS_DOUBLE) ||
       !SDDS_DefineSimpleParameter(&SDDSout, "S34Sigma", "m", SDDS_DOUBLE) ||
       !SDDS_DefineSimpleParameter(&SDDSout, "S44Sigma", "", SDDS_DOUBLE) ||
       !SDDS_DefineSimpleParameter(&SDDSout, "S13Sigma", "m$a2$n", SDDS_DOUBLE) ||
       !SDDS_DefineSimpleParameter(&SDDSout, "S14Sigma", "m", SDDS_DOUBLE) ||
       !SDDS_DefineSimpleParameter(&SDDSout, "S23Sigma", "m", SDDS_DOUBLE) ||
       !SDDS_DefineSimpleParameter(&SDDSout, "S24Sigma", "", SDDS_DOUBLE) ||
       !SDDS_DefineSimpleParameter(&SDDSout, "S16Sigma", "m", SDDS_DOUBLE) ||
       !SDDS_DefineSimpleParameter(&SDDSout, "S26Sigma", "", SDDS_DOUBLE) ||
       !SDDS_DefineSimpleParameter(&SDDSout, "S36Sigma", "m", SDDS_DOUBLE) ||
       !SDDS_DefineSimpleParameter(&SDDSout, "S46Sigma", "", SDDS_DOUBLE) ||
       !SDDS_DefineSimpleParameter(&SDDSout, "goodFits", "", SDDS_LONG) ||
       !SDDS_DefineSimpleParameter(&SDDSout, "averageFitPoints", "", SDDS_DOUBLE) ||
       !SDDS_DefineSimpleParameter(&SDDSout, "exSigma", "m", SDDS_DOUBLE) ||
       !SDDS_DefineSimpleParameter(&SDDSout, "eySigma", "m", SDDS_DOUBLE) ||
       !SDDS_DefineSimpleParameter(&SDDSout, "SdeltaSigma", "", SDDS_DOUBLE) ||
       !SDDS_DefineSimpleParameter(&SDDSout, "betaxSigma", "m", SDDS_DOUBLE) ||
       !SDDS_DefineSimpleParameter(&SDDSout, "betaySigma", "m", SDDS_DOUBLE) ||
       !SDDS_DefineSimpleParameter(&SDDSout, "alphaxSigma", "", SDDS_DOUBLE) ||
       !SDDS_DefineSimpleParameter(&SDDSout, "alphaySigma", "", SDDS_DOUBLE) ||
       !SDDS_DefineSimpleParameter(&SDDSout, "etaxSigma", "m", SDDS_DOUBLE) ||
       !SDDS_DefineSimpleParameter(&SDDSout, "etaySigma", "m", SDDS_DOUBLE) ||
       !SDDS_DefineSimpleParameter(&SDDSout, "etaxpSigma", "", SDDS_DOUBLE) ||
       !SDDS_DefineSimpleParameter(&SDDSout, "etaypSigma", "", SDDS_DOUBLE))) {
    SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors | SDDS_VERBOSE_PrintErrors);
  }

  if (SDDS_GetColumnIndex(&SDDSin, "S11") < 0 ||
      SDDS_GetColumnIndex(&SDDSin, "S33") < 0 ||
      SDDS_GetColumnIndex(&SDDSin, "S13") < 0 ||
      SDDS_GetColumnIndex(&SDDSin, "R11") < 0 ||
      SDDS_GetColumnIndex(&SDDSin, "R12") < 0 ||
      SDDS_GetColumnIndex(&SDDSin, "R33") < 0 ||
      SDDS_GetColumnIndex(&SDDSin, "R34") < 0 ||
      SDDS_GetColumnIndex(&SDDSin, "R16") < 0 ||
      SDDS_GetColumnIndex(&SDDSin, "R26") < 0)
    SDDS_Bomb("input file missing required quantities. Need S11, S33, S13, R11, R12, R33, R34, R16, and R26.");

  if (nErrorSets > 2 &&
      (SDDS_GetColumnIndex(&SDDSin, "S11StDev") < 0 ||
       SDDS_GetColumnIndex(&SDDSin, "S33StDev") < 0 ||
       SDDS_GetColumnIndex(&SDDSin, "S13StDev") < 0))
    SDDS_Bomb("input file missing required quantities. Need S11StDev, S33StDev, and S13StDev.");

  if (!SDDS_WriteLayout(&SDDSout))
    SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors | SDDS_VERBOSE_PrintErrors);

  /* suppress compiler warnings */
  R11 = R12 = R33 = R34 = R16 = R26 = NULL;
  S11 = S33 = S13 = NULL;
  S11StDev = S33StDev = S13StDev = NULL;

  while (SDDS_ReadTable(&SDDSin) > 0) {
    nConfigs = SDDS_CountRowsOfInterest(&SDDSin);
    if (nConfigs < 4)
      continue;

    if (!(R11 = SDDS_GetColumn(&SDDSin, "R11")) ||
        !(R12 = SDDS_GetColumn(&SDDSin, "R12")) ||
        !(R33 = SDDS_GetColumn(&SDDSin, "R33")) ||
        !(R34 = SDDS_GetColumn(&SDDSin, "R34")) ||
        !(R16 = SDDS_GetColumn(&SDDSin, "R16")) ||
        !(R26 = SDDS_GetColumn(&SDDSin, "R26")))
      SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors | SDDS_VERBOSE_PrintErrors);

    if (!(S11 = SDDS_GetColumn(&SDDSin, "S11")) ||
        !(S33 = SDDS_GetColumn(&SDDSin, "S33")) ||
        !(S13 = SDDS_GetColumn(&SDDSin, "S13")))
      SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors | SDDS_VERBOSE_PrintErrors);

    if (nErrorSets > 2 &&
        (!(S11StDev = SDDS_GetColumn(&SDDSin, "S11StDev")) ||
         !(S33StDev = SDDS_GetColumn(&SDDSin, "S33StDev")) ||
         !(S13StDev = SDDS_GetColumn(&SDDSin, "S13StDev"))))
      SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors | SDDS_VERBOSE_PrintErrors);

    if (!SDDS_StartPage(&SDDSout, nConfigs) || !SDDS_CopyColumns(&SDDSout, &SDDSin))
      SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors | SDDS_VERBOSE_PrintErrors);

    if (x_scale != 1 || y_scale != 1) {
      long i;
      for (i = 0; i < nConfigs; i++) {
        S11[i] *= sqr(x_scale);
        S33[i] *= sqr(y_scale);
        S13[i] *= x_scale * y_scale;
        if (S11StDev)
          S11StDev[i] *= sqr(x_scale);
        if (S33StDev)
          S33StDev[i] *= sqr(y_scale);
        if (S13StDev)
          S13StDev[i] *= x_scale * y_scale;
      }
    }

    for (i = 0; i < 5; i++)
      for (j = 0; j < 5; j++)
        SijSum[i][j] = SijSum2[i][j] = 0;
    for (i = 0; i < 4; i++)
      for (j = 0; j < 4; j++)
        SbetaSum[i][j] = SbetaSum2[i][j] = 0;
    for (i = 0; i < 2; i++)
      emitSum[i] = emitSum2[i] = betaSum[i] = betaSum2[i] = alphaSum[i] = alphaSum2[i] =
        etaSum[i] = etapSum[i] = 0;

    S11e = calloc(nConfigs, sizeof(*S11e));
    S13e = calloc(nConfigs, sizeof(*S13e));
    S33e = calloc(nConfigs, sizeof(*S33e));
    S11Fit = calloc(nConfigs, sizeof(*S11Fit));
    S13Fit = calloc(nConfigs, sizeof(*S13Fit));
    S33Fit = calloc(nConfigs, sizeof(*S33Fit));
    S11FitSum = calloc(nConfigs, sizeof(*S11FitSum));
    S13FitSum = calloc(nConfigs, sizeof(*S13FitSum));
    S33FitSum = calloc(nConfigs, sizeof(*S33FitSum));
    S11FitSum2 = calloc(nConfigs, sizeof(*S11FitSum2));
    S13FitSum2 = calloc(nConfigs, sizeof(*S13FitSum2));
    S33FitSum2 = calloc(nConfigs, sizeof(*S33FitSum2));

    for (iSet = 0; iSet < nErrorSets; iSet++) {
      double beta[2], alpha[2], eta[2], etap[2], emit[2];
      double Sbeta[4][4];
      short goodResult, trialLimit;

      for (i = 0; i < 4; i++)
        for (j = 0; j < 4; j++)
          Sbeta[i][j] = 0;

      goodResult = 0;
      if (nErrorSets > 2)
        trialLimit = 100;
      else
        trialLimit = 1;
      do {
        for (iConfig = 0; iConfig < nConfigs; iConfig++) {
          if (nErrorSets > 2) {
            S11e[iConfig] = S11[iConfig] + gauss_rn_lim(0.0, S11StDev[iConfig], 2, random_1) - sqr(x_resol);
            S13e[iConfig] = S13[iConfig] + gauss_rn_lim(0.0, S13StDev[iConfig], 2, random_1) - SIGN(S13[iConfig]) * x_resol * y_resol;
            S33e[iConfig] = S33[iConfig] + gauss_rn_lim(0.0, S33StDev[iConfig], 2, random_1) - sqr(y_resol);
          } else {
            S11e[iConfig] = S11[iConfig] - sqr(x_resol);
            S13e[iConfig] = S13[iConfig] - SIGN(S13[iConfig]) * x_resol * y_resol;
            S33e[iConfig] = S33[iConfig] - sqr(y_resol);
          }
        }

        for (i = 0; i < 5; i++)
          for (j = 0; j < 5; j++)
            Sij[i][j] = -1;

        solveForSigmaMatrix(1, 1, S11e, nErrorSets > 2 ? S11StDev : NULL, R11, R12, R33, R34, R16, R26, nConfigs, Sij, S11Fit, energySpread);
        solveForSigmaMatrix(3, 3, S33e, nErrorSets > 2 ? S33StDev : NULL, R11, R12, R33, R34, R16, R26, nConfigs, Sij, S33Fit, energySpread);
        solveForSigmaMatrix(1, 3, S13e, nErrorSets > 2 ? S13StDev : NULL, R11, R12, R33, R34, R16, R26, nConfigs, Sij, S13Fit, energySpread);
        goodResult = 1;
        for (i = 0; i < 5; i++) {
          if (Sij[i][i] < 0) {
            goodResult = 0;
            break;
          }
          for (j = 0; j < 5; j++) {
            if (isnan(Sij[i][j]) || isinf(Sij[i][j])) {
              goodResult = 0;
              break;
            }
          }
        }
        trialLimit--;
      } while (trialLimit > 0 && !goodResult);
      if (!goodResult)
        SDDS_Bomb("Failed to find any good solutions.");

      for (i = 0; i < nConfigs; i++) {
        S11FitSum[i] += S11Fit[i];
        S13FitSum[i] += S13Fit[i];
        S33FitSum[i] += S33Fit[i];
        S11FitSum2[i] += sqr(S11Fit[i]);
        S13FitSum2[i] += sqr(S13Fit[i]);
        S33FitSum2[i] += sqr(S33Fit[i]);
      }

      for (i = 0; i < 2; i++) {
        eta[i] = Sij[2 * i + 0][4] / Sij[4][4];
        etap[i] = Sij[2 * i + 1][4] / Sij[4][4];
        etaSum[i] += eta[i];
        etaSum2[i] += sqr(eta[i]);
        etapSum[i] += etap[i];
        etapSum2[i] += sqr(etap[i]);

        Sbeta[2 * i + 0][2 * i + 0] = Sij[2 * i + 0][2 * i + 0] - sqr(eta[i]) * Sij[4][4];
        Sbeta[2 * i + 1][2 * i + 1] = Sij[2 * i + 1][2 * i + 1] - sqr(etap[i]) * Sij[4][4];
        Sbeta[2 * i + 0][2 * i + 1] = Sij[2 * i + 0][2 * i + 1] - eta[i] * etap[i] * Sij[4][4];
        Sbeta[2 * i + 1][2 * i + 0] = Sbeta[2 * i + 0][2 * i + 1];

        emit[i] = sqrt(Sbeta[2 * i + 0][2 * i + 0] * Sbeta[2 * i + 1][2 * i + 1] - sqr(Sbeta[2 * i + 0][2 * i + 1]));
        beta[i] = Sbeta[2 * i + 0][2 * i + 0] / emit[i];
        alpha[i] = -Sbeta[2 * i + 0][2 * i + 1] / emit[i];
      }

      for (i = 0; i < 5; i++)
        for (j = 0; j < 5; j++) {
          SijSum[i][j] += Sij[i][j];
          SijSum2[i][j] += sqr(Sij[i][j]);
        }
      for (i = 0; i < 4; i++)
        for (j = 0; j < 4; j++) {
          SbetaSum[i][j] += Sbeta[i][j];
          SbetaSum2[i][j] += sqr(Sbeta[i][j]);
        }
      for (i = 0; i < 2; i++) {
        emitSum[i] += emit[i];
        emitSum2[i] += sqr(emit[i]);
        betaSum[i] += beta[i];
        betaSum2[i] += sqr(beta[i]);
        alphaSum[i] += alpha[i];
        alphaSum2[i] += sqr(alpha[i]);
      }
    }

    if (!SDDS_StartPage(&SDDSout, nConfigs))
      SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors | SDDS_VERBOSE_PrintErrors);

    if (!SDDS_SetParameters(&SDDSout, SDDS_SET_BY_NAME | SDDS_PASS_BY_VALUE,
                            "ex", emitSum[0] / nErrorSets,
                            "ey", emitSum[1] / nErrorSets,
                            "betax", betaSum[0] / nErrorSets,
                            "betay", betaSum[1] / nErrorSets,
                            "alphax", alphaSum[0] / nErrorSets,
                            "alphay", alphaSum[1] / nErrorSets,
                            "etax", etaSum[0] / nErrorSets,
                            "etay", etaSum[1] / nErrorSets,
                            "etaxp", etapSum[0] / nErrorSets,
                            "etayp", etapSum[1] / nErrorSets,
                            "Sdelta", sqrt(SijSum[4][4] / nErrorSets),
                            NULL))
      SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors | SDDS_VERBOSE_PrintErrors);

    if (nErrorSets > 2 &&
        !SDDS_SetParameters(&SDDSout, SDDS_SET_BY_NAME | SDDS_PASS_BY_VALUE,
                            "exSigma", sqrt(emitSum2[0] / nErrorSets - sqr(emitSum[0] / nErrorSets)),
                            "eySigma", sqrt(emitSum2[1] / nErrorSets - sqr(emitSum[1] / nErrorSets)),
                            "betaxSigma", sqrt(betaSum2[0] / nErrorSets - sqr(betaSum[0] / nErrorSets)),
                            "betaySigma", sqrt(betaSum2[1] / nErrorSets - sqr(betaSum[1] / nErrorSets)),
                            "alphaxSigma", sqrt(alphaSum2[0] / nErrorSets - sqr(alphaSum[0] / nErrorSets)),
                            "alphaySigma", sqrt(alphaSum2[1] / nErrorSets - sqr(alphaSum[1] / nErrorSets)),
                            "etaxSigma", sqrt(etaSum2[0] / nErrorSets - sqr(etaSum[0] / nErrorSets)),
                            "etaySigma", sqrt(etaSum2[1] / nErrorSets - sqr(etaSum[1] / nErrorSets)),
                            "etaxpSigma", sqrt(etapSum2[0] / nErrorSets - sqr(etapSum[0] / nErrorSets)),
                            "etaypSigma", sqrt(etapSum2[1] / nErrorSets - sqr(etapSum[1] / nErrorSets)),
                            "SdeltaSigma", sqrt(sqrt(SijSum2[4][4] / nErrorSets - sqr(SijSum[4][4] / nErrorSets))),
                            NULL))
      SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors | SDDS_VERBOSE_PrintErrors);

    if (!SDDS_SetParameters(&SDDSout, SDDS_SET_BY_NAME | SDDS_PASS_BY_VALUE,
                            "S11", SijSum[0][0] / nErrorSets, "S12", SijSum[0][1] / nErrorSets, "S22", SijSum[1][1] / nErrorSets,
                            NULL))
      SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors | SDDS_VERBOSE_PrintErrors);

    if (nErrorSets > 2 &&
        !SDDS_SetParameters(&SDDSout, SDDS_SET_BY_NAME | SDDS_PASS_BY_VALUE,
                            "S11Sigma", sqrt(SijSum2[0][0] / nErrorSets - sqr(SijSum[0][0] / nErrorSets)),
                            "S12Sigma", sqrt(SijSum2[0][1] / nErrorSets - sqr(SijSum[0][1] / nErrorSets)),
                            "S22Sigma", sqrt(SijSum2[1][1] / nErrorSets - sqr(SijSum[1][1] / nErrorSets)),
                            NULL))
      SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors | SDDS_VERBOSE_PrintErrors);

    if (!SDDS_SetParameters(&SDDSout, SDDS_SET_BY_NAME | SDDS_PASS_BY_VALUE,
                            "S16", SijSum[0][4] / nErrorSets, "S26", SijSum[1][4] / nErrorSets, "S66", SijSum[4][4] / nErrorSets,
                            NULL))
      SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors | SDDS_VERBOSE_PrintErrors);

    if (nErrorSets > 2 &&
        !SDDS_SetParameters(&SDDSout, SDDS_SET_BY_NAME | SDDS_PASS_BY_VALUE,
                            "S16Sigma", sqrt(SijSum2[0][4] / nErrorSets - sqr(SijSum[0][4] / nErrorSets)),
                            "S26Sigma", sqrt(SijSum2[1][4] / nErrorSets - sqr(SijSum[0][4] / nErrorSets)),
                            "S66Sigma", sqrt(SijSum2[4][4] / nErrorSets - sqr(SijSum[4][4] / nErrorSets)),
                            NULL))
      SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors | SDDS_VERBOSE_PrintErrors);

    if (!SDDS_SetParameters(&SDDSout, SDDS_SET_BY_NAME | SDDS_PASS_BY_VALUE,
                            "S33", SijSum[2][2] / nErrorSets, "S34", SijSum[2][3] / nErrorSets, "S44", SijSum[3][3] / nErrorSets,
                            NULL))
      SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors | SDDS_VERBOSE_PrintErrors);

    if (nErrorSets > 2 &&
        !SDDS_SetParameters(&SDDSout, SDDS_SET_BY_NAME | SDDS_PASS_BY_VALUE,
                            "S33Sigma", sqrt(SijSum2[2][2] / nErrorSets - sqr(SijSum[2][2] / nErrorSets)),
                            "S34Sigma", sqrt(SijSum2[2][3] / nErrorSets - sqr(SijSum[2][3] / nErrorSets)),
                            "S44Sigma", sqrt(SijSum2[3][3] / nErrorSets - sqr(SijSum[3][3] / nErrorSets)),
                            NULL))
      SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors | SDDS_VERBOSE_PrintErrors);

    if (!SDDS_SetParameters(&SDDSout, SDDS_SET_BY_NAME | SDDS_PASS_BY_VALUE,
                            "S13", SijSum[0][2] / nErrorSets, "S14", SijSum[0][3] / nErrorSets,
                            "S23", SijSum[1][2] / nErrorSets, "S24", SijSum[1][3] / nErrorSets,
                            NULL))
      SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors | SDDS_VERBOSE_PrintErrors);

    if (nErrorSets > 2 &&
        !SDDS_SetParameters(&SDDSout, SDDS_SET_BY_NAME | SDDS_PASS_BY_VALUE,
                            "S13Sigma", sqrt(SijSum2[0][2] / nErrorSets - sqr(SijSum[0][2] / nErrorSets)),
                            "S14Sigma", sqrt(SijSum2[0][3] / nErrorSets - sqr(SijSum[0][3] / nErrorSets)),
                            "S23Sigma", sqrt(SijSum2[1][2] / nErrorSets - sqr(SijSum[1][2] / nErrorSets)),
                            "S24Sigma", sqrt(SijSum2[1][3] / nErrorSets - sqr(SijSum[1][3] / nErrorSets)),
                            NULL))
      SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors | SDDS_VERBOSE_PrintErrors);

    if (!SDDS_SetParameters(&SDDSout, SDDS_SET_BY_NAME | SDDS_PASS_BY_VALUE,
                            "S36", SijSum[2][4] / nErrorSets, "S46", SijSum[3][4] / nErrorSets,
                            NULL))
      SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors | SDDS_VERBOSE_PrintErrors);

    if (nErrorSets > 2 &&
        !SDDS_SetParameters(&SDDSout, SDDS_SET_BY_NAME | SDDS_PASS_BY_VALUE,
                            "S36Sigma", sqrt(SijSum2[2][4] / nErrorSets - sqr(SijSum[2][4] / nErrorSets)),
                            "S46Sigma", sqrt(SijSum2[3][4] / nErrorSets - sqr(SijSum[3][4] / nErrorSets)),
                            NULL))
      SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors | SDDS_VERBOSE_PrintErrors);

    for (i = 0; i < nConfigs; i++) {
      S11FitSum[i] /= nErrorSets;
      S13FitSum[i] /= nErrorSets;
      S33FitSum[i] /= nErrorSets;
    }
    if (!SDDS_SetColumn(&SDDSout, SDDS_SET_BY_NAME, S11, nConfigs, "S11Data") ||
        !SDDS_SetColumn(&SDDSout, SDDS_SET_BY_NAME, S33, nConfigs, "S33Data") ||
        !SDDS_SetColumn(&SDDSout, SDDS_SET_BY_NAME, S13, nConfigs, "S13Data") ||
        !SDDS_SetColumn(&SDDSout, SDDS_SET_BY_NAME, S11FitSum, nConfigs, "S11Fit") ||
        !SDDS_SetColumn(&SDDSout, SDDS_SET_BY_NAME, S33FitSum, nConfigs, "S33Fit") ||
        !SDDS_SetColumn(&SDDSout, SDDS_SET_BY_NAME, S13FitSum, nConfigs, "S13Fit") ||
        !SDDS_WritePage(&SDDSout))
      SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors | SDDS_VERBOSE_PrintErrors);

    free(R11);
    free(R12);
    free(R33);
    free(R34);
    free(R16);
    free(R26);
    free(S11);
    free(S33);
    free(S13);
    if (nErrorSets > 2) {
      free(S11StDev);
      free(S33StDev);
      free(S13StDev);
    }
    free(S11Fit);
    free(S33Fit);
    free(S13Fit);
    free(S11FitSum);
    free(S33FitSum);
    free(S13FitSum);
    free(S11FitSum2);
    free(S33FitSum2);
    free(S13FitSum2);
    free(S11e);
    free(S13e);
    free(S33e);
  }

  if (!SDDS_Terminate(&SDDSin) || !SDDS_Terminate(&SDDSout))
    SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors | SDDS_VERBOSE_PrintErrors);

  return (0);
}

void solveForSigmaMatrix(int i, int j,
                         double *SMeasured, /* measured Sij with error added */
                         double *SStDev,    /* standard deviation in values */
                         double *R11, double *R12, double *R33, double *R34, double *R16, double *R26,
                         long nConfigs,
                         double Sij[5][5], double *SFit,
                         double energySpread) {
  /* Sets up problem to solve M = P.F for one of the three cases described in the comment at the top
     of this file.
  */
  long iConfig, nUnknowns;
  MATRIX *P, *F, *M, *K, *FSigma;

  if (i == j && i == 1) {
    nUnknowns = 6;
    if (energySpread > 0)
      nUnknowns = 5;
  } else if (i == j && i == 3) {
    nUnknowns = 3;
  } else if (i == 1 && j == 3) {
    nUnknowns = 6;
  } else {
    bomb("Invalid parameters for solveForSigmaMatrix. Seek expert help.", NULL);
    return;
  }

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
        P->a[iConfig][3] = 2 * R11[iConfig] * R16[iConfig];
        P->a[iConfig][4] = 2 * R12[iConfig] * R16[iConfig];
        if (energySpread <= 0)
          P->a[iConfig][5] = sqr(R16[iConfig]);
      } else {
        P->a[iConfig][0] = sqr(R33[iConfig]);
        P->a[iConfig][1] = 2 * R33[iConfig] * R34[iConfig];
        P->a[iConfig][2] = sqr(R34[iConfig]);
      }
    } else {
      P->a[iConfig][0] = R11[iConfig] * R33[iConfig];
      P->a[iConfig][1] = R11[iConfig] * R34[iConfig];
      P->a[iConfig][2] = R12[iConfig] * R33[iConfig];
      P->a[iConfig][3] = R12[iConfig] * R34[iConfig];
      P->a[iConfig][4] = R16[iConfig] * R33[iConfig];
      P->a[iConfig][5] = R16[iConfig] * R34[iConfig];
    }
  }

  /* Setup covariance matrix */
  m_zero(K);
  if (SStDev) {
    for (iConfig = 0; iConfig < nConfigs; iConfig++)
      K->a[iConfig][iConfig] = 1 / sqr(2 * SStDev[iConfig]);
  } else {
    for (iConfig = 0; iConfig < nConfigs; iConfig++)
      K->a[iConfig][iConfig] = 1;
  }

  for (iConfig = 0; iConfig < nConfigs; iConfig++) {
    M->a[iConfig][0] = SMeasured[iConfig];
    if (energySpread > 0 && i == 1 && j == 1)
      M->a[iConfig][0] -= sqr(R16[iConfig]) * sqr(energySpread);
  }

  solve_normal_form(F, FSigma, P, M, K, SFit);
  if (i == j) {
    if (i == 1) {
      Sij[0][0] = F->a[0][0];
      Sij[0][1] = F->a[1][0];
      Sij[1][1] = F->a[2][0];
      Sij[0][4] = F->a[3][0];
      Sij[1][4] = F->a[4][0];
      if (energySpread <= 0)
        Sij[4][4] = F->a[5][0];
      else
        Sij[4][4] = sqr(energySpread);
    } else {
      Sij[2][2] = F->a[0][0];
      Sij[2][3] = F->a[1][0];
      Sij[3][3] = F->a[2][0];
    }
  } else {
    Sij[0][2] = F->a[0][0];
    Sij[0][3] = F->a[1][0];
    Sij[1][2] = F->a[2][0];
    Sij[1][3] = F->a[3][0];
    Sij[2][4] = F->a[4][0];
    Sij[3][4] = F->a[5][0];
  }

  m_free(&P);
  m_free(&F);
  m_free(&FSigma);
  m_free(&M);
  m_free(&K);

  return;
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

void set_up_covariance_matrix(MATRIX *K, double *sigma, double *uncert, long nConfigs, long equal_weights)
/* actually, K is the inverse of the covariance matrix */
{
  long i_config;

  for (i_config = 0; i_config < nConfigs; i_config++) {
    if (sigma[i_config] == 0 || uncert[i_config] == 0 || equal_weights)
      K->a[i_config][i_config] = 1;
    else
      K->a[i_config][i_config] = 1. / sqr(2 * sigma[i_config] * uncert[i_config]);
  }
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
