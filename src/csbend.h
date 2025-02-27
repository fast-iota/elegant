#ifndef _CSBEND_H_
#define _CSBEND_H_

#define EXSQRT(value, order) (order==0?sqrt(value):(1+0.5*((value)-1)))

#ifdef __cplusplus
extern "C" {
#endif

/* instances of these globals are in csbend.c */
/*
extern long negativeWarningsLeft;
extern long dipoleFringeWarning;
*/

// variables can get inlined better in csbend if not needed for gpu
// will be declared static in csbend.c
#if defined(HAVE_GPU)
extern long expansionOrder1;  /* order of expansion+1 */
extern long hasSkew, hasNormal;
extern double rho0, rho_actual, rad_coef, isrConstant;
extern double meanPhotonsPerRadian0, meanPhotonsPerMeter0, normalizedCriticalEnergy0;
extern long distributionBasedRadiation, includeOpeningAngle;
extern long photonCount;
extern double energyCount, radiansTotal;
extern double **Fx_xy, **Fy_xy;
extern short refTrajectoryMode;
extern long refTrajectoryPoints;
extern double **refTrajectoryData;
#endif

#if !defined(PARALLEL)
/* to avoid problems with HP parallel compiler */
extern unsigned long multipoleKicksDone ;
#endif

#define DERBENEV_CRITERION_DISABLE 0
#define DERBENEV_CRITERION_EVAL 1
#define DERBENEV_CRITERION_ENFORCE 2
#define N_DERBENEV_CRITERION_OPTIONS 3

typedef struct {
  unsigned long lastMode;
#define CSRDRIFT_STUPAKOV          0x0001UL
#define CSRDRIFT_SALDIN54          0x0002UL
#define CSRDRIFT_OVERTAKINGLENGTH  0x0004UL
#define CSRDRIFT_ATTENUATIONLENGTH 0x0008UL
#define CSRDRIFT_SPREAD            0x0010UL
  long bins, valid;
  double dctBin, s0, ds0, zLast, z0;
  double S11, S12, S22;
  double *dGamma;
  double rho, bendingAngle, Po, perc68BunchLength, perc90BunchLength, peakToPeakWavelength, rmsBunchLength;
  /* for Saldin eq 54 (NIM A 398 (1997) 373-394) mode: */
  FILE *fpSaldin;
  long nSaldin;
  double *FdNorm;   /* Saldin's Fd(sh, x)/Fd(sh, 0), sh = bunch-length*gamma^3/rho */
  double *xSaldin;  /* distance from end of bend */
  double lastFdNorm; /* last value obtained from interpolation */
  /* for Stupakov mode */
  long SGOrder, SGHalfWidth, SGDerivHalfWidth, SGDerivOrder;
  double binRangeFactor;
  double GSConstant, MPCharge;
  char *StupakovOutput;
  SDDS_DATASET SDDS_Stupakov;
  long StupakovFileActive, StupakovOutputInterval;
  long trapazoidIntegration;
  double lowFrequencyCutoff0, lowFrequencyCutoff1;
  double highFrequencyCutoff0, highFrequencyCutoff1;
  long clipNegativeBins;
  long wffValues;
  double *wffFreqValue, *wffRealFactor, *wffImagFactor;
} CSR_LAST_WAKE;

extern CSR_LAST_WAKE csrWake;

#ifdef __cplusplus
}
#endif

#endif
