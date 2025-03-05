/*************************************************************************\
* Copyright (c) 2002 The University of Chicago, as Operator of Argonne
* National Laboratory.
* Copyright (c) 2002 The Regents of the University of California, as
* Operator of Los Alamos National Laboratory.
* This file is distributed subject to a Software License Agreement found
* in the file LICENSE that is included with this distribution. 
\*************************************************************************/

#include "mdb.h"
#include "SDDS.h"
#include "scan.h"

#define SET_PIPE 0
#define SET_DE1 1
#define SET_COLUMNS 2
#define N_OPTIONS 3

#define NPPSIGMA 6

char *option[N_OPTIONS] = {"pipe", "de1", "columns"};

char *USAGE = "sddsecon [<input>] [<output>] [-pipe=[in][,out]]\n\
               -columns=<energy-name>,<spectrum-name>,<output-spectrum-name>\n\
               -de1=<value>\n\n\
<input>        The input file should contain a column with photon energy (ev,keV)\n\
               and a spectrum column.\n\
-de1           Detector energy resolution at lower energy (ev,keV)\n\
RESTRICTIONS:\n\
 The energy points need to be equally spaced in ascending order.\n\
 Due to the convolution performed, some data points will be set to\n\
 zero at each end of the spectrum. The energy range with zeroed data\n\
 will be bigger at the end of the spectrum due to the transformation used,\n\
 and it depends on the value of the beam energy spread. The larger the\n\
 beam energy spread, the larger the range with zeroed data.\n\
PROCEDURE:\n\
 The spectrum is convolved with a Gaussian function with a constant width.\n\n\
Program by Robert Soliday. ("__DATE__")\n";

double* econ(double *e, double *spec, long rows, double de1);
long checkMonotonicity(double *indepValue, long rows);

int main(int argc, char **argv)
{
  SCANNED_ARG *s_arg;
  long i_arg;

  double de1=0;
  char *input=NULL, *output=NULL, *eUnits=NULL, *specUnits=NULL;
  char *eName=NULL, *specName=NULL, *specOutputName=NULL;
  unsigned long pipeFlags=0;
  long noWarnings=0, tmpfile_used=0, rows;

  SDDS_DATASET SDDS_input, SDDS_output;
  double *e, *spec, *specResults=NULL;

  SDDS_RegisterProgramName(argv[0]);
  argc = scanargs(&s_arg, argc, argv);
  if (argc<4) {
    fprintf(stderr, "%s", USAGE);
    return(1);
  }
  /* parse the command line */
  for (i_arg=1; i_arg<argc; i_arg++) {
    if (s_arg[i_arg].arg_type==OPTION) {
      switch (match_string(s_arg[i_arg].list[0], option, N_OPTIONS, 0)) {
      case SET_PIPE:
        if (!processPipeOption(s_arg[i_arg].list+1, s_arg[i_arg].n_items-1, &pipeFlags))
          SDDS_Bomb("invalid -pipe syntax");
        break;
      case SET_COLUMNS:
        if (s_arg[i_arg].n_items!=4)
          SDDS_Bomb("invalid -columns syntax");
        eName = s_arg[i_arg].list[1];
        specName = s_arg[i_arg].list[2];
        specOutputName = s_arg[i_arg].list[3];
        break;
      case SET_DE1:
	if (s_arg[i_arg].n_items!=2) {
	  fprintf(stderr, "error: invalid -de1 syntax\n");
	  return(1);
	}
	if (sscanf(s_arg[i_arg].list[1], "%lf", &de1) != 1) {
	  fprintf(stderr, "error: invalid -de1 syntax or value\n");
	  return(1);
	}
	break;
      default:
	fprintf(stderr, "error: unknown switch: %s\n", s_arg[i_arg].list[0]);
	return(1);
      }
    } else {
      if (input == NULL) {
	input = s_arg[i_arg].list[0];
      } else if (output == NULL) {
	output = s_arg[i_arg].list[0];
      } else {
	fprintf(stderr, "too many filenames\n");
	return(1);
      }
    }
  }
  if (de1 == 0) {
    fprintf(stderr, "error: invalid -de1 syntax or value\n");
    return(1);
  }
  if ((eName == NULL) || (specName == NULL)) {
    fprintf(stderr, "error: invalid -columns syntax or value\n");
    return(1);
  }
  processFilenames("sddsecon", &input, &output, pipeFlags, noWarnings, &tmpfile_used);

  /* open the input file */
  if (!SDDS_InitializeInput(&SDDS_input, input)) {
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
    return(1);
  }
  if (!SDDS_GetColumnInformation(&SDDS_input, "units", &eUnits, SDDS_GET_BY_NAME, eName)) {
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
    return(1);
  }
  if (!SDDS_GetColumnInformation(&SDDS_input, "units", &specUnits, SDDS_GET_BY_NAME, specName)) {
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
    return(1);
  }
  if (SDDS_ReadPage(&SDDS_input) != 1) {
    fprintf(stderr, "error: no data found in input file\n");
    return(1);
  }
  rows = SDDS_RowCount(&SDDS_input);

  if (!(e = SDDS_GetColumnInDoubles(&SDDS_input, eName)) ||
      !(spec = SDDS_GetColumnInDoubles(&SDDS_input, specName))) {
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
    return(1);
  }
  /* close the input file */
  if (!SDDS_Terminate(&SDDS_input)) {
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
    return(1);
  }

  specResults = econ(e, spec, rows, de1);
  if (specResults == NULL) {
    return(1);
  }
  /* open the output file */
  if (!SDDS_InitializeOutput(&SDDS_output, SDDS_BINARY, 0, NULL, NULL, output)) {
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
    return(1);
  }
  if (!SDDS_DefineSimpleColumn(&SDDS_output, eName, eUnits, SDDS_DOUBLE)) {
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
    return(1);
  }
  if (!SDDS_DefineSimpleColumn(&SDDS_output, specName, specUnits, SDDS_DOUBLE)) {
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
    return(1);
  }
  if (!SDDS_DefineSimpleColumn(&SDDS_output, specOutputName, specUnits, SDDS_DOUBLE)) {
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
    return(1);
  }
  if (!SDDS_WriteLayout(&SDDS_output)) {
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
    return(1);
  }
  if (!SDDS_StartTable(&SDDS_output, rows)) {
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
    return(1);
  }
  if (!SDDS_SetColumn(&SDDS_output, SDDS_SET_BY_NAME, e, rows, eName)) {
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
    return(1);
  }
  if (!SDDS_SetColumn(&SDDS_output, SDDS_SET_BY_NAME, spec, rows, specName)) {
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
    return(1);
  }
  if (!SDDS_SetColumn(&SDDS_output, SDDS_SET_BY_NAME, specResults, rows, specOutputName)) {
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
    return(1);
  }
  if (!SDDS_WriteTable(&SDDS_output)) {
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
    return(1);
  }
  /* close the input file */
  if (!SDDS_Terminate(&SDDS_output)) {
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
    return(1);
  }

  free(e);
  free(spec);
  free(specResults);
  free_scanargs(&s_arg,argc);
  return(0);
}

double* econ(double *e, double *spec, long rows, double de1) {
  double *cf, *cg, *specg;
  long ng, i, monotonicity;

  double q, sigf;
  unsigned long interpCode=0;
  OUTRANGE_CONTROL aboveRange, belowRange;

  long nsigma, nsigma2, n_pts, t, k2;
  double conv, xvar, *gaus, gausSum, *specc, *specResults;

  aboveRange.flags = belowRange.flags = OUTRANGE_ABORT;

  cf = malloc(sizeof(double) * rows);
  q = (rows - 1) / log(e[rows-1] / e[0]);
  for (i = 0; i < rows; i++) {
    cf[i] = q * log(e[i]/e[0]); /* new distribution 'cf' with unequal step size */
  }
  sigf = 2.0 * de1 * q; 
  if (sigf < .01) {
    /* Return original spec array because the value of sigf is small */
    specResults = malloc(sizeof(double) * rows);
    for (i = 0; i < rows; i++) {
      specResults[i] = spec[i];
    }
    return(specResults);
  }
  /* Number of channels for new distribution with equidistant step size */
  ng = ceil(NPPSIGMA * rows / sigf);
  if (rows > ng)
    ng = rows;
  /* Generate channel distribution with equidistant step size */
  cg = malloc(sizeof(double) * ng);
  specg = malloc(sizeof(double) * ng);
  for (i = 0; i < ng; i++) {
    cg[i] = i * (rows - 1) / (ng - 1);
  }

  /* New function interpolated into new distribution with equidistant step sizes */
  if ((monotonicity=checkMonotonicity(cf, rows))==0) {
    fprintf(stderr, "independent data values do not change monotonically or repeated independent values exist\n");
    return(NULL);
  }
  for (i = 0; i < ng; i++) {
    specg[i] = interpolate(spec, cf, rows, cg[i], &belowRange, &aboveRange, 1, &interpCode, monotonicity);
    if (interpCode) {
      if (interpCode&OUTRANGE_ABORT) {
        fprintf(stderr, "error: can't do interpolation, value %lf is out of range.\n", cg[i]);
        return(NULL);
      }
    }
  }

  /* Prep for convolution with a gaussian with constant width */
  nsigma = 3.0; /* nsigma is a factor that determines how many points to use for the gaussian kernel. 
                   the exact number of points used depends on the value of sigf and the range of X values
                   as well, but will be roughly equal to 2*NSIGMA*sigf */
  nsigma2 = nsigma * 2;
  conv = (cg[ng-1] - cg[0]) / (ng - 1); /* conversion, units/point */
  n_pts = ceil(nsigma2 * sigf / conv); /* number of points */
  /* Restrict n_pts */
  if (2 > n_pts)
    n_pts = 2;
  if ((ng - 2) < n_pts)
    n_pts = ng - 2;
  if ((long)(n_pts/(long)2 * 2) == n_pts)
    n_pts += 1;
  /* Create gaussian */
  gaus = malloc(sizeof(double) * n_pts);
  gausSum = 0;
  for (i = 0; i < n_pts; i++) {
    xvar = (i / (n_pts - 1.0) - .5) * n_pts;
    gaus[i] = exp(-.5 * pow((xvar / (sigf / conv)),2));
    gausSum += gaus[i];
  }
  /* Convolve with gaussian whose width is sigf.
     http://www.exelisvis.com/docs/CONVOL.html */
  specc = malloc(sizeof(double) * ng);
  k2 = n_pts / (long)2;
  for (t=0; t < ng; t++) {
    specc[t]=0;
    for (i=0; i <= n_pts - 1; i++) {
      if ((k2 <= t) && (t <= ng - k2 - 1)) {
        specc[t] += specg[t + i - k2] * gaus[i];
      }
    }
    specc[t] = specc[t] / gausSum;
  }

  /* Revert the cf distribution */
  specResults = malloc(sizeof(double) * rows);
  if ((monotonicity=checkMonotonicity(cg, ng))==0) {
    fprintf(stderr, "independent data values do not change monotonically or repeated independent values exist\n");
    return(NULL);
  }
  for (i = 0; i < rows; i++) {
    specResults[i] = interpolate(specc, cg, ng, cf[i], &belowRange, &aboveRange, 1, &interpCode, monotonicity);
    if (interpCode) {
      if (interpCode&OUTRANGE_ABORT) {
        fprintf(stderr, "error: can't do interpolation, value %lf is out of range.\n", cf[i]);
        return(NULL);
      }
    }
  }

  free(cf);
  free(cg);
  free(specg);
  free(gaus);
  free(specc);
  return(specResults);
}

long checkMonotonicity(double *indepValue, long rows)
{
  if (rows==1)
    return 1;
  if (indepValue[rows-1]>indepValue[0]) {
    while (--rows>0)
      if (indepValue[rows]<=indepValue[rows-1])
        return 0;
    return 1;
  }
  else {
    while (--rows>0)
      if (indepValue[rows]>=indepValue[rows-1])
        return 0;
    return -1;
  }
}
