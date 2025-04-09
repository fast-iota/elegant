#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

void brighte(int i, double alpha, double cosphi, double sinphi,
             double *s0, double *s1, double *s2, double *s3, double kx, double ky);

void hunt(const double *X, int N, double XP, int *jlo);

#define E_SZ 40001
#define A_SZ 400
#define B_SZ (A_SZ / 4 + 1)
#define ZERO (0.0)
#define ONE (1.0)
#define TWO (2.0)
#define C_ANG_M (1.0e-10)
#define C_M2_MM2 (1.0e6)
#define CV (1.0 / (8.0 * PI * PI) * C_ANG_M * C_M2_MM2)
#define IDWMAX 128

#define C_ 2.99792458e8      /* Speed of light [m/s]         */
static const double MEE = 0.51099906;       /* Electron rest mass [MeV]     */
#define EC 1.60217733e-19    /* Elementary charge [C]        */
#define H_ 6.6260755e-34     /* Planck~s constant [J s]      */
#define HBAR 1.05457266e-34  /* h/2Pi [J s]                  */
#define EPSZ 8.854187817e-12 /* Permittivity of vacuum [F/m] */
#define PI 3.14159265358979323846
static const double PIHALF = 1.57079632679489661923;
static const double TWOPI = 6.28318530717958647692;

/* Conversion factors */
static const double C_EVANG = (H_ * C_ / EC * 1.0e10);
static const double C_M_MM = 1.0e3;
static const double C_CM_ANG = 1.0e8;
static const double C_MRAD_RAD = 1.0e-3;
static const double C_MA_A = 1.0e-3;
static const double C_CM_M = 1.0e-2;

static const double FINE_STRUCTURE_CONST =
  (EC * 1.0e19 * EC * 1.0e19) / (4.0 * PI * EPSZ * 1.0e38 * HBAR * C_);
static const double BW = 1.0e-3; /* 0.1% bandwidth */
static const double EPS = 1.0e-6;
static const double EPSE = 1.0e-8;
static const double EPSK = 1.0e-4;

typedef struct {
  int METHOD, IHARM, N, NSIGMA;
  int NALPHA;
  int NE, NE1, NE2, NEU;
  double KX, KY, K3, NPI, PE, GAMMA, ER, EW, DEW, LEN;
  double FAC, C1;
  double AP2MIN, AP2MAX, AP2CNT, ARGMAX, FU, FV, SIGX2, SIGY2;
  double CALPHA2, DE;
  int I1[E_SZ], I2[E_SZ];
  double E[E_SZ], EU[E_SZ], SPEC0[E_SZ], SPEC1[E_SZ], SPEC2[E_SZ], SPEC3[E_SZ];
  int NW;
  double HW[E_SZ], DW;
  double HE[E_SZ], HA[E_SZ];
  int INDEX_PHI[A_SZ], NPHI4, NPHI, NPHI_BRIGHT;
  double COSPHI[A_SZ], SINPHI[A_SZ], S2SIGN[A_SZ], DPHI;
} GLOBAL_VALS;

/**
 * @brief Computes the squared sinc function of the given argument.
 *
 * The function calculates sinc(arg) = sin(arg)/arg, squared.
 * For very small arguments (|arg| < 1.0e-6), it returns 1.0 to avoid numerical instability.
 *
 * @param arg The input value for which to compute the sinc function.
 * @return The squared sinc of the input argument.
 */
static double sinc(double arg) {
  double x = fabs(arg);
  if (x > 1.0e-6) {
    double sl = sin(arg) / arg;
    return sl * sl;
  }
  return 1.0;
}

/**
 * @brief Performs the energy convolution step for the spectral distribution.
 *
 * This function convolutes the energy distribution by applying a smoothing kernel to the input spectral array \c SP1.
 * The convolution is performed between the indices \c gv->NE1 - 1 and \c gv->NE2 - 1.
 *
 * @param SP1 Pointer to the input/output spectral array to be convoluted.
 * @param gv Pointer to a \c GLOBAL_VALS structure containing global parameters and arrays.
 */
void convolute_energy_estep(double *SP1, GLOBAL_VALS *gv) {
  double SP2[E_SZ];
  int IE, IW;

  for (IE = gv->NE1 - 1; IE < gv->NE2; IE++) {
    SP2[IE] = 0.0;
    for (IW = 0; IW < gv->NW; IW++) {
      int idx = (IE + 1) + gv->NE1 - (IW + 1) - 1;
      if (idx >= 0 && idx < E_SZ) {
        SP2[IE] += SP1[idx] * gv->HW[IW];
      }
    }
  }
  for (IE = gv->NE1 - 1; IE < gv->NE2; IE++) {
    SP1[IE] = SP2[IE] * gv->DW;
  }
}

/**
 * @brief Performs the energy convolution step for velocity space in the spectral distribution.
 *
 * This function convolutes the spectral distribution \c SP1 by integrating over energy using a velocity step method.
 * It updates \c SP1 based on the convolution with the hardware weightings and other global parameters.
 *
 * @param SP1 Pointer to the input/output spectral array to be convoluted.
 * @param gv Pointer to a \c GLOBAL_VALS structure containing global parameters and arrays.
 */
void convolute_energy_vstep(double *SP1, GLOBAL_VALS *gv) {
  double SP2[E_SZ];
  int IEU, J1 = 0, J2 = 0;

  for (IEU = 0; IEU < gv->NEU; IEU++)
    SP2[IEU] = 0.0;

  for (IEU = 0; IEU < gv->NEU; IEU++) {
    double EUP = gv->EU[IEU];
    double EU1 = EUP - gv->EW;
    double EU2 = EUP + gv->EW;
    hunt(gv->E, gv->NE, EU1, &J1);
    hunt(gv->E, gv->NE, EU2, &J2);

    if (J2 == -1 || J1 == gv->NE - 1) {
      SP2[IEU] = 0.0;
    } else if (J1 == -1) {
      if (J2 < gv->NE - 1) {
        double H2 = EU2 - gv->E[J2];
        double P = H2 / gv->HE[J2];
        double S2 = (ONE - P) * SP1[J2] + P * SP1[J2 + 1];
        double ARG = gv->PE * (EUP - EU2);
        double SUM = S2 * sinc(ARG) * H2 / TWO;
        ARG = gv->PE * (EUP - gv->E[J2]);
        SUM += SP1[J2] * sinc(ARG) * (H2 + gv->HE[J2 - 1] / TWO);
        for (int i = 0; i < J2; i++) {
          ARG = gv->PE * (EUP - gv->E[i]);
          SUM += SP1[i] * sinc(ARG) * gv->HA[i];
        }
        SP2[IEU] = SUM / gv->DEW;
      } else {
        double SUM = 0.0;
        for (int i = 0; i < gv->NE; i++) {
          double ARG = gv->PE * (EUP - gv->E[i]);
          SUM += SP1[i] * sinc(ARG) * gv->HA[i];
        }
        SP2[IEU] = SUM / gv->DEW;
      }
    } else if (J2 < gv->NE - 1) {
      if (J1 == J2) {
        SP2[IEU] = 0.0;
      } else {
        double H1 = gv->E[J1 + 1] - EU1;
        double P = H1 / gv->HE[J1];
        double S1 = P * SP1[J1] + (ONE - P) * SP1[J1 + 1];
        double ARG = gv->PE * (EUP - EU1);
        double SUM = S1 * sinc(ARG) * H1 / TWO;
        ARG = gv->PE * (EUP - gv->E[J1 + 1]);
        SUM += SP1[J1 + 1] * sinc(ARG) * (H1 + gv->HE[J1 + 1]) / TWO;

        double H2 = EU2 - gv->E[J2];
        P = H2 / gv->HE[J2];
        double S2 = (ONE - P) * SP1[J2] + P * SP1[J2 + 1];
        ARG = gv->PE * (EUP - EU2);
        SUM += S2 * sinc(ARG) * H2 / TWO;
        ARG = gv->PE * (EUP - gv->E[J2]);
        SUM += SP1[J2] * sinc(ARG) * (H2 + gv->HE[J2 - 1] / TWO);

        for (int i = J1 + 2; i < J2; i++) {
          ARG = gv->PE * (EUP - gv->E[i]);
          SUM += SP1[i] * sinc(ARG) * gv->HA[i];
        }
        SP2[IEU] = SUM / gv->DEW;
      }
    } else {
      double H1 = gv->E[J1 + 1] - EU1;
      double P = H1 / gv->HE[J1];
      double S1 = P * SP1[J1] + (ONE - P) * SP1[J1 + 1];
      double ARG = gv->PE * (EUP - EU1);
      double SUM = S1 * sinc(ARG) * H1 / TWO;
      ARG = gv->PE * (EUP - gv->E[J1 + 1]);
      SUM += SP1[J1 + 1] * sinc(ARG) * (H1 + gv->HE[J1 + 1]) / TWO;
      for (int i = J1 + 2; i < gv->NE; i++) {
        ARG = gv->PE * (EUP - gv->E[i]);
        SUM += SP1[i] * sinc(ARG) * gv->HA[i];
      }
      SP2[IEU] = SUM / gv->DEW;
    }
  }
  for (IEU = 0; IEU < gv->NEU; IEU++) {
    SP1[IEU] = SP2[IEU];
  }
}

/**
 * @brief Convolves the distribution over alpha and theta angles.
 *
 * This function performs the convolution of the brightness arrays BR0_, BR1_, BR2_, and BR3_ with the alpha and theta parameters.
 * It updates the accumulated results RA0_, RA1_, RA2_, RA3_ and the count ICOUNT based on the calculation method.
 *
 * @param ICALC Calculation mode indicator (e.g., finite-N or infinite-N).
 * @param CONST_ Constant scaling factor used in the convolution.
 * @param ALPHA Pointer to an array of alpha values.
 * @param THETA Pointer to an array of theta values.
 * @param DALPHA Step size for alpha.
 * @param EPS_ Epsilon value for convergence criteria.
 * @param BR0_ Two-dimensional array representing the first brightness component.
 * @param BR1_ Two-dimensional array representing the second brightness component.
 * @param BR2_ Two-dimensional array representing the third brightness component.
 * @param BR3_ Two-dimensional array representing the fourth brightness component.
 * @param RA0_ Pointer to accumulate the first result component.
 * @param RA1_ Pointer to accumulate the second result component.
 * @param RA2_ Pointer to accumulate the third result component.
 * @param RA3_ Pointer to accumulate the fourth result component.
 * @param ICOUNT Pointer to an integer tracking the number of counts or iterations.
 * @param gv Pointer to a \c GLOBAL_VALS structure containing global parameters and arrays.
 */
void convolute_distribution(int ICALC, double CONST_, double *ALPHA, double *THETA,
                            double DALPHA, double EPS_, double BR0_[B_SZ][B_SZ],
                            double BR1_[B_SZ][B_SZ], double BR2_[B_SZ][B_SZ],
                            double BR3_[B_SZ][B_SZ], double *RA0_, double *RA1_,
                            double *RA2_, double *RA3_, int *ICOUNT, GLOBAL_VALS *gv) {
  double SL = CONST_ * DALPHA * gv->DPHI;
  int IC, ID;
  if (ICALC == 1) { /* Finite-N */
    double SUM0 = 0.0, SUM1 = 0.0, SUM2 = 0.0, SUM3 = 0.0;
    for (IC = 0; IC < gv->NALPHA; IC++) {
      double DELTA0 = 0.0, DELTA1 = 0.0, DELTA2 = 0.0, DELTA3 = 0.0;
      for (ID = 0; ID < gv->NPHI4; ID++) {
        double U = THETA[IC] * gv->COSPHI[ID];
        double V = THETA[IC] * gv->SINPHI[ID];
        double ARG = U * U * gv->FU + V * V * gv->FV;
        if (ARG < gv->ARGMAX) {
          double P = exp(-ARG);
          DELTA0 += P * BR0_[gv->INDEX_PHI[ID]][IC];
          DELTA1 += P * BR1_[gv->INDEX_PHI[ID]][IC];
          DELTA2 += P * BR2_[gv->INDEX_PHI[ID]][IC] * gv->S2SIGN[ID];
          DELTA3 += P * BR3_[gv->INDEX_PHI[ID]][IC];
        }
      }
      SUM0 += DELTA0 * ALPHA[IC];
      SUM1 += DELTA1 * ALPHA[IC];
      SUM2 += DELTA2 * ALPHA[IC];
      SUM3 += DELTA3 * ALPHA[IC];
    }
    SUM0 *= SL;
    SUM1 *= SL;
    SUM2 *= SL;
    SUM3 *= SL;
    (*RA0_) += SUM0;
    (*RA1_) += SUM1;
    (*RA2_) += SUM2;
    (*RA3_) += SUM3;
    if (SUM0 > EPS_ * (*RA0_))
      *ICOUNT = 0;
  } else { /* Infinite-N */
    double DELTA0 = 0.0, DELTA1 = 0.0, DELTA2 = 0.0, DELTA3 = 0.0;
    for (ID = 0; ID < gv->NPHI4; ID++) {
      double U = THETA[0] * gv->COSPHI[ID];
      double V = THETA[0] * gv->SINPHI[ID];
      double ARG = U * U * gv->FU + V * V * gv->FV;
      if (ARG < gv->ARGMAX) {
        double P = exp(-ARG);
        DELTA0 += P * BR0_[gv->INDEX_PHI[ID]][0];
        DELTA1 += P * BR1_[gv->INDEX_PHI[ID]][0];
        DELTA2 += P * BR2_[gv->INDEX_PHI[ID]][0] * gv->S2SIGN[ID];
        DELTA3 += P * BR3_[gv->INDEX_PHI[ID]][0];
      }
    }
    DELTA0 *= SL;
    DELTA1 *= SL;
    DELTA2 *= SL;
    DELTA3 *= SL;
    (*RA0_) += DELTA0;
    (*RA1_) += DELTA1;
    (*RA2_) += DELTA2;
    (*RA3_) += DELTA3;
    if (DELTA0 > EPS_ * (*RA0_))
      *ICOUNT = 0;
  }
}

/**
 * @brief Computes the brightness arrays based on calculation mode and alpha parameters.
 *
 * This function calculates the brightness components BR0, BR1, BR2, and BR3 for each angle ID and alpha index IC.
 * It applies a sinc modulation based on the argument derived from the alpha values.
 *
 * @param ICALC Calculation mode indicator (e.g., finite-N or infinite-N).
 * @param I_ Current harmonic index.
 * @param R_ Scaling factor related to energy.
 * @param ALPHAI Current alpha value.
 * @param ALPHA2I Current alpha squared value.
 * @param ALPHA_ Pointer to an array of alpha values.
 * @param BR0 Two-dimensional array to store the first brightness component.
 * @param BR1 Two-dimensional array to store the second brightness component.
 * @param BR2 Two-dimensional array to store the third brightness component.
 * @param BR3 Two-dimensional array to store the fourth brightness component.
 * @param gv Pointer to a \c GLOBAL_VALS structure containing global parameters and arrays.
 */
void brightness_array(int ICALC, int I_, double R_, double ALPHAI, double ALPHA2I,
                      double *ALPHA_, double BR0[B_SZ][B_SZ], double BR1[B_SZ][B_SZ],
                      double BR2[B_SZ][B_SZ], double BR3[B_SZ][B_SZ],
                      GLOBAL_VALS *gv) {
  int IC, ID;
  if (ICALC == 1) { /* Finite-N */
    for (IC = 0; IC < gv->NALPHA; IC++) {
      double AL2 = ALPHA_[IC] * ALPHA_[IC];
      double ARG = gv->NPI * (AL2 - ALPHA2I) / R_;
      double H = sinc(ARG);
      for (ID = 0; ID < gv->NPHI_BRIGHT; ID++) {
        brighte(I_, ALPHA_[IC], gv->COSPHI[ID], gv->SINPHI[ID],
                &(BR0[ID][IC]), &(BR1[ID][IC]), &(BR2[ID][IC]), &(BR3[ID][IC]), gv->KX, gv->KY);
        BR0[ID][IC] *= H;
        BR1[ID][IC] *= H;
        BR2[ID][IC] *= H;
        BR3[ID][IC] *= H;
      }
    }
  } else { /* Infinite-N */
    for (ID = 0; ID < gv->NPHI_BRIGHT; ID++) {
      brighte(I_, ALPHAI, gv->COSPHI[ID], gv->SINPHI[ID],
              &(BR0[ID][0]), &(BR1[ID][0]), &(BR2[ID][0]), &(BR3[ID][0]), gv->KX, gv->KY);
    }
  }
}

/**
 * @brief Computes the spectral distribution based on global parameters.
 *
 * This function calculates the spectral distribution for each energy level, considering harmonics and applying convolution steps.
 * It populates the SPEC0, SPEC1, SPEC2, and SPEC3 arrays within the \c GLOBAL_VALS structure.
 *
 * @param IERROR Pointer to an integer flag for error reporting (0 = no error, -1 = error).
 * @param gv Pointer to a \c GLOBAL_VALS structure containing global parameters and arrays.
 */
void spectral_distribution(int32_t *IERROR, GLOBAL_VALS *gv) {
  *IERROR = 0;
  int LE1 = 0, LE2 = 0;

  for (int IE = 0; IE < gv->NE; IE++) {
    double BR0[B_SZ][B_SZ], BR1[B_SZ][B_SZ], BR2[B_SZ][B_SZ], BR3[B_SZ][B_SZ];
    // Initialize spectral and index arrays
    gv->SPEC0[IE] = gv->SPEC1[IE] = gv->SPEC2[IE] = gv->SPEC3[IE] = 0.0;
    gv->I1[IE] = gv->I2[IE] = 0;

    double R_ = gv->ER / gv->E[IE];
    double DA2 = (gv->METHOD == 1) ? (gv->CALPHA2 * R_ / gv->N) : 0.0;

    // Calculate harmonic range
    int IMIN = (int)((gv->AP2MIN + gv->K3 - ((gv->METHOD == 1) ? DA2 : 0.0)) / R_ + 1.0);
    int IMAX = (int)((gv->AP2MAX + gv->K3 + ((gv->METHOD == 1) ? DA2 : 0.0)) / R_);

    if (IMAX < IMIN)
      continue; // No valid harmonics for this energy

    LE1 = 1;

    // Check if specific harmonic is within range
    if (gv->IHARM > 0 && (gv->IHARM < IMIN || gv->IHARM > IMAX))
      continue;

    LE2 = 1;

    // Initialize accumulation variables
    double RA0 = 0.0, RA1 = 0.0, RA2 = 0.0, RA3 = 0.0;

    // Determine harmonic indices based on IHARM
    int IH1, IH2;
    if (gv->IHARM > 0) {
      IH1 = IH2 = gv->IHARM;
    } else if (gv->IHARM < 0) {
      IH1 = IH2 = IMIN;
    } else {
      IH1 = IMIN;
      IH2 = IMAX;
    }

    gv->I1[IE] = IH1;
    gv->I2[IE] = IH2;
    int ICOUNT = 0;

    // Loop over harmonics
    for (int I = IH1; I <= IH2; I++) {
      ICOUNT++;
      double ALPHA2I = R_ * I - gv->K3;

      double ALPHAI, DEL, DALPHA;
      double CONST_;
      double ALPHArr[A_SZ], THETAarr[A_SZ];

      if (gv->METHOD == 1) {
        // METHOD 1: Handle multiple alpha values
        double AL2m = fmax(ALPHA2I - DA2, 0.0);
        double ALm = sqrt(AL2m);
        double ALM = sqrt(ALPHA2I + DA2);
        DALPHA = (ALM - ALm) / gv->NALPHA;
        double SL = ALm + DALPHA / 2.0;

        for (int IC = 0; IC < gv->NALPHA; IC++) {
          ALPHArr[IC] = SL + (IC * DALPHA);
          THETAarr[IC] = ALPHArr[IC] / gv->GAMMA;
        }

        ALPHAI = (ALPHA2I > 0.0) ? sqrt(ALPHA2I) : 0.0;
        DEL = (ALPHA2I < 0.0) ? 0.0 : ALPHA2I * gv->N / R_;
      } else {
        // Other METHODS: Single alpha value
        ALPHAI = sqrt(ALPHA2I);
        ALPHArr[0] = ALPHAI;
        THETAarr[0] = ALPHAI / gv->GAMMA;

        DALPHA = R_ / (2 * gv->N);
        DEL = (ALPHA2I < 0.0) ? 0.0 : ALPHA2I * gv->N / R_;
      }

      // Calculate SIGR2 based on DEL
      double tmp, SIGR2;
      if (DEL < 2.15) {
        // tmp = 1.29 + 1.229 * pow(DEL - 0.8, 2);
        tmp = 1.29 + 1.229 * (DEL - 0.8) * (DEL - 0.8);
        SIGR2 = tmp * tmp;
      } else {
        SIGR2 = 5.81 * DEL;
      }
      SIGR2 *= CV * C_EVANG / gv->E[IE] * gv->LEN;

      // Calculate CONST_
      CONST_ = gv->C1 / (TWOPI * sqrt((gv->SIGX2 + SIGR2) * (gv->SIGY2 + SIGR2)));

      // Populate brightness arrays
      brightness_array(gv->METHOD, I, R_, ALPHAI, ALPHA2I, ALPHArr, BR0, BR1, BR2, BR3, gv);

      // Perform convolution
      convolute_distribution(gv->METHOD, CONST_, ALPHArr, THETAarr, DALPHA,
                             EPS, BR0, BR1, BR2, BR3,
                             &RA0, &RA1, &RA2, &RA3, &ICOUNT, gv);

      // Exit loop if multiple counts are detected
      if (ICOUNT > 1) {
        gv->I2[IE] = I;
        break;
      }
    }

    // Accumulate spectral results
    gv->SPEC0[IE] = gv->FAC * RA0;
    gv->SPEC1[IE] = gv->FAC * RA1;
    gv->SPEC2[IE] = 0.0; // Symmetry
    gv->SPEC3[IE] = gv->FAC * RA3;
  }

  // Error handling based on LE1 and LE2 flags
  if (!LE1 && !LE2) {
    printf("&SPECTRAL_DISTRIBUTION-E-HARMERR, no harmonics reachable\n");
    *IERROR = -1;
    return;
  }
  if (LE1 && !LE2) {
    printf("&SPECTRAL_DISTRIBUTION-E-HARMERR, harmonic not in range\n");
    *IERROR = -1;
    return;
  }

  // Apply final convolution based on METHOD
  if (gv->METHOD == 4) {
    convolute_energy_vstep(gv->SPEC0, gv);
    convolute_energy_vstep(gv->SPEC1, gv);
    convolute_energy_vstep(gv->SPEC2, gv);
    convolute_energy_vstep(gv->SPEC3, gv);
  } else if (gv->METHOD == 14) {
    convolute_energy_estep(gv->SPEC0, gv);
    convolute_energy_estep(gv->SPEC1, gv);
    convolute_energy_estep(gv->SPEC2, gv);
    convolute_energy_estep(gv->SPEC3, gv);
  }
}

/**
 * @brief Computes the on-axis brilliance spectral distribution for a given energy and beam parameters.
 *
 * This function initializes global parameters based on input beam characteristics, sets up energy and angle arrays,
 * and calls the spectral distribution computation. It returns the energy array and corresponding on-axis brilliance
 * as output arrays. Error flags are set in case of invalid input or computation issues.
 *
 * @param ENERGY Electron beam energy [eV].
 * @param CUR Beam current [A].
 * @param SIGX Horizontal beam size [m].
 * @param SIGY Vertical beam size [m].
 * @param SIGX1 Horizontal beam divergence [rad].
 * @param SIGY1 Vertical beam divergence [rad].
 * @param PERIOD Undulator period [m].
 * @param N_F Number of points in the energy grid.
 * @param KX_F Horizontal wave number component.
 * @param KY_F Vertical wave number component.
 * @param EMINU Minimum energy to consider [eV].
 * @param EMAXU Maximum energy to consider [eV].
 * @param NEU_F Number of energy points in the output.
 * @param METHOD_F Calculation method identifier:
 *        - 0: Dejus' infinite-N + convolution
 *        - 1: Dejus' method (finite-N)
 *        - 2: Walker's infinite-N
 *        - 3: Walker's finite-N
 * @param E_F Output array for energy values [eV].
 * @param SPEC0_F Output array for on-axis brilliance.
 * @param NE_F Pointer to an integer to store the number of points returned.
 * @param IERROR Pointer to an integer error flag (0 = success, -1 = error).
 */
void usb(
         double ENERGY, double CUR, double SIGX, double SIGY, double SIGX1, double SIGY1,
         double PERIOD, int N_F, double KX_F, double KY_F,
         double EMINU, double EMAXU, int NEU_F, int METHOD_F,
         double *E_F,     /* output array (energy) */
         double *SPEC0_F, /* output array (on-axis brilliance) */
         int32_t *NE_F,   /* number of points returned */
         int32_t *IERROR  /* error flag; 0=ok, -1=error */
         ) {
  double COMEGA, D2, G2, LAMDAR, K2, E1Z, SIGU2, SIGV2, SIGU, SIGV, XE, YE, E1MIN, E1MAX, EMIN, EMAX, ARG, SL;
  int I, IE1, IE2, IE, IP1, IP2, IP, ISIGN, ID;
  GLOBAL_VALS gv;
  *IERROR = 0;

  /* -------------------------------------------------------------------
     Copy the Fortran logic for assigning N, KX, KY, NEU, etc.
     ------------------------------------------------------------------- */
  gv.N = N_F;
  gv.KX = KX_F;
  gv.KY = KY_F;
  gv.NEU = NEU_F;

  /* Map METHOD_F to internal METHOD codes (mirroring Fortran code) */
  if (METHOD_F == 0)
    gv.METHOD = 4; /* Dejus~ infinite-N + convolution */
  else if (METHOD_F == 1)
    gv.METHOD = 4;
  else if (METHOD_F == 2)
    gv.METHOD = 14; /* Walker~s infinite-N */
  else if (METHOD_F == 3)
    gv.METHOD = 1; /* Walker~s finite-N */
  else
    goto Error930;

  if (fabs(gv.KX - gv.KY) >= EPSK && (gv.KX >= EPSK && gv.KY >= EPSK)) {
    /* ~Kx, Ky must match or one must be zero.~ */
    goto Error902;
  }

  /* -------------------------------------------------------------------
     Default values (from Fortran)
     ------------------------------------------------------------------- */
  gv.NPHI = 20;
  gv.NALPHA = 15;
  gv.NSIGMA = 3;
  gv.IHARM = 0; /* all harmonics included => IHARM=0 */
  gv.CALPHA2 = 2.0;
  COMEGA = 8.0;
  /* In the Fortran code, D=1 and used in D^2, but effectively just 1.0^2 */
  D2 = 1.0;

  /* -------------------------------------------------------------------
     Calculate ~global~ parameters
     ------------------------------------------------------------------- */
  {
    gv.GAMMA = (ENERGY / MEE) * 1.0e3;
    G2 = gv.GAMMA * gv.GAMMA;
    LAMDAR = PERIOD * C_CM_ANG / (2.0 * G2); /* [Angstroms]   */
    gv.ER = C_EVANG / LAMDAR;                /* [eV]          */
    K2 = gv.KX * gv.KX + gv.KY * gv.KY;
    gv.K3 = 1.0 + K2 / 2.0;
    E1Z = gv.ER / gv.K3;             /* 1st harmonic energy [eV]            */
    gv.LEN = gv.N * PERIOD * C_CM_M; /* device length [m] */
    gv.NPI = gv.N * PI;
  }

  /* -------------------------------------------------------------------
     Beam emittance
     ------------------------------------------------------------------- */
  gv.SIGX2 = SIGX * SIGX;
  gv.SIGY2 = SIGY * SIGY;
  SIGU2 = (SIGX1 * SIGX1) * (C_MRAD_RAD * C_MRAD_RAD);
  SIGV2 = (SIGY1 * SIGY1) * (C_MRAD_RAD * C_MRAD_RAD);
  SIGU = sqrt(SIGU2);
  SIGV = sqrt(SIGV2);
  gv.FU = (SIGU2 != 0.0) ? (0.5 / SIGU2) : 0.0;
  gv.FV = (SIGV2 != 0.0) ? (0.5 / SIGV2) : 0.0;

  /* -------------------------------------------------------------------
     Determine min, max emission angles
     ------------------------------------------------------------------- */
  gv.AP2MIN = 0.0;
  gv.AP2CNT = 0.0;
  XE = gv.NSIGMA * SIGU;
  YE = gv.NSIGMA * SIGV;
  gv.AP2MAX = G2 * (XE * XE + YE * YE);
  gv.ARGMAX = (gv.NSIGMA * gv.NSIGMA) / 2.0;

  /* -------------------------------------------------------------------
     Walker / Dejus approach for energy scale
     ------------------------------------------------------------------- */
  if (gv.METHOD == 4) {                /* Dejus~ method */
    gv.GAMMA = (ENERGY / MEE) * 1.0e3; /* re-compute if needed */
    G2 = gv.GAMMA * gv.GAMMA;
    /* E1Z was computed above. */
    E1MIN = E1Z * gv.K3 / (gv.K3 + gv.AP2MAX);
    E1MAX = E1Z * gv.K3 / (gv.K3 + gv.AP2MIN);
    gv.DEW = (E1Z / N_F) * gv.K3 / (gv.K3 + gv.AP2CNT);
    gv.EW = COMEGA * gv.DEW;
    gv.PE = PI / gv.DEW;

    /* Extend lower/upper limit */
    EMIN = EMINU - gv.EW;
    EMAX = EMAXU + gv.EW;
    if (EMIN < EPS)
      EMIN = EPS;

    /* Find IE1, IE2 in the same style as Fortran, etc. */
    I = 1;
    {
      double EP_local = I * E1MAX;
      while (EP_local < EMIN) {
        I++;
        EP_local = I * E1MAX;
      }
      IE1 = I;

      if ((IE1 * E1MIN) > EMIN) {
        EMIN = IE1 * E1MIN;
      }

      EP_local = IE1 * E1MIN;
      while (EP_local < EMAX) {
        I++;
        EP_local = I * E1MIN;
      }
      IE2 = I - 1;

      if (IE2 < IE1)
        goto Error910; /* EMAX <= EMIN => no harmonics */
    }

    /* Variable-step energy mesh creation */
    {
      int IEcount = 0;
      gv.E[IEcount] = EMIN;
      I = IE1;
      {
        double EP_local = I * E1MAX;
        while (gv.E[IEcount] < EMAX) {
          int IDW_local = 1;
          double EF_local = gv.DEW / IDW_local;
          double DW_local = 2.0 * gv.EW / (double)(64); /* NOMEGA=64 default? */

          /* reduce step if needed */
          while (gv.E[IEcount] > (EP_local - EF_local)) {
            IDW_local *= 2;
            if (IDW_local > IDWMAX)
              EF_local = 0.0;
            else {
              DW_local = (2.0 * gv.EW) / (double)(64 * IDW_local);
              EF_local = gv.DEW / (double)(IDW_local);
            }
          }

          while (gv.E[IEcount] < (EP_local - EPSE) && (gv.E[IEcount] < (EMAX - EPSE))) {
            if (gv.E[IEcount] > (EP_local - EF_local)) {
              IDW_local *= 2;
              if (IDW_local > IDWMAX)
                EF_local = 0.0;
              else {
                DW_local = (2.0 * gv.EW) / (double)(64 * IDW_local);
                EF_local = gv.DEW / (double)(IDW_local);
              }
            }
            IEcount++;
            if (IEcount >= E_SZ)
              goto Error900;
            gv.HE[IEcount - 1] = DW_local;
            if (IEcount == 1)
              gv.HA[IEcount - 1] = gv.HE[IEcount - 1] / 2.0;
            else
              gv.HA[IEcount - 1] = (gv.HE[IEcount - 2] + gv.HE[IEcount - 1]) / 2.0;

            gv.E[IEcount] = gv.E[IEcount - 1] + gv.HE[IEcount - 1];
          }

          /* Force a point near EP (within EPSE) */
          SL = (EP_local < EMAX) ? EP_local : EMAX;
          gv.E[IEcount] = SL - EPSE;
          gv.HE[IEcount - 1] = gv.E[IEcount] - gv.E[IEcount - 1];
          gv.HA[IEcount - 1] = (gv.HE[IEcount - 2] + gv.HE[IEcount - 1]) / 2.0;
          IEcount++;
          if (IEcount >= E_SZ)
            goto Error900;
          gv.HE[IEcount - 1] = 2.0 * EPSE;
          gv.HA[IEcount - 1] = (gv.HE[IEcount - 2] + gv.HE[IEcount - 1]) / 2.0;
          gv.E[IEcount] = gv.E[IEcount - 1] + gv.HE[IEcount - 1];

          /* Next harmonic index */
          I++;
          {
            double DEF_local = I * E1MIN - gv.E[IEcount];
            if (DEF_local > 0.0) {
              IEcount++;
              if (IEcount >= E_SZ)
                goto Error900;
              gv.HE[IEcount - 1] = DEF_local;
              gv.HA[IEcount - 1] = (gv.HE[IEcount - 2] + gv.HE[IEcount - 1]) / 2.0;
              gv.E[IEcount] = gv.E[IEcount - 1] + gv.HE[IEcount - 1];
            }
          }

          EP_local = I * E1MAX;
        }
        gv.NE = IEcount;
      }

      gv.HE[gv.NE] = 0.0;
      gv.HA[gv.NE] = gv.HE[gv.NE - 1] / 2.0;
      if (gv.NEU == 0)
        gv.NEU = (int)(gv.NE / 100.0 + 1.0) * 100; /* default in Fortran */

      gv.DE = (EMAXU - EMINU) / (double)gv.NEU; /* user~s desired interval */
      gv.NEU = gv.NEU + 1;

      for (IE = 0; IE < gv.NEU; IE++)
        gv.EU[IE] = EMINU + (IE)*gv.DE;
    }
    gv.NE1 = 1;
    gv.NE2 = gv.NE;
  } else if (gv.METHOD == 14) { /* Walker~s infinite-N approach */
    gv.GAMMA = (ENERGY / MEE) * 1.0e3;
    G2 = gv.GAMMA * gv.GAMMA; /* re-ensure it */
    /* E1Z computed above, K3 above, etc. */
    gv.DEW = (E1Z / N_F) * gv.K3 / (gv.K3 + gv.AP2CNT);
    gv.EW = COMEGA * gv.DEW;
    if (gv.NEU > 0) {
      /* re-define NOMEGA from Fortran approach: (EMAXU-EMINU)/NEU, etc. */
      SL = (EMAXU - EMINU) / (double)gv.NEU;
      int NOMEGA_local = (int)(2.0 * gv.EW / SL + 1.0);
      NOMEGA_local = (NOMEGA_local / 2) * 2; /* force even */
      if (NOMEGA_local < 16) {
        NOMEGA_local = 16;
        printf("\n&USB-W-ISMALL, NOMEGA less than 16 => reset to 16\n");
      }
      gv.NW = NOMEGA_local + 1;
    } else {
      gv.NW = 64 + 1; /* default if NEU=0, say 64 from Fortran */
    }

    gv.DW = (2.0 * COMEGA) / (double)(gv.NW - 1);
    for (I = 0; I < gv.NW; I++) {
      ARG = (-COMEGA + (I)*gv.DW) * PI;
      if (fabs(ARG) > EPS) {
        SL = sin(ARG) / ARG;
        gv.HW[I] = SL * SL;
      } else {
        gv.HW[I] = 1.0;
      }
    }

    /* Step in energy, set up E[] array */
    gv.DE = gv.DW * gv.DEW;
    EMIN = EMINU - gv.EW;
    if (EMIN < EPS)
      EMIN = EPS;
    EMAX = EMAXU + gv.EW;
    gv.NE = (int)((EMAX - EMIN) / gv.DE) + 1;
    gv.NE = gv.NE + 1;            /* # points in E array */
    gv.NE1 = (gv.NW - 1) / 2 + 1; /* from Fortran: NOMEGA/2+1 => center index */
    gv.NE2 = gv.NE + 1 - gv.NE1;

    if (gv.NE > E_SZ)
      goto Error905;
    for (IE = 0; IE < gv.NE; IE++)
      gv.E[IE] = EMIN + (IE)*gv.DE;

  } else { /* METHOD == 1 => Walker~s finite-N method (or other fallback) */
    gv.NE = gv.NEU;
    if (gv.NE > 0) {
      gv.DE = (EMAXU - EMINU) / (double)gv.NE;
    } else {
      fprintf(stderr, "Error: METHOD==1 and neks value must be greater than 0\n");
      exit(1);
    }
    gv.NE = gv.NE + 1;
    gv.NE1 = 1;
    gv.NE2 = gv.NE + 1 - gv.NE1;
    if (gv.NE > E_SZ)
      goto Error905;

    for (IE = 0; IE < gv.NE; IE++)
      gv.E[IE] = EMINU + (IE)*gv.DE;
  }

  /* -------------------------------------------------------------------
     Set up arrays for angles
     ------------------------------------------------------------------- */
  gv.NPHI = 20;
  gv.NPHI4 = gv.NPHI;
  gv.NPHI_BRIGHT = gv.NPHI;

  gv.DPHI = PIHALF / gv.NPHI;
  ISIGN = +1;
  SL = ISIGN * gv.DPHI / TWO;
  gv.INDEX_PHI[0] = 0;
  gv.S2SIGN[0] = +ONE;
  gv.COSPHI[0] = cos(SL);
  gv.SINPHI[0] = sin(SL);
  for (ID = 1; ID < gv.NPHI4; ID++) {
    ARG = SL + (ID * gv.DPHI);
    gv.INDEX_PHI[ID] = gv.INDEX_PHI[ID - 1] + ISIGN;
    gv.COSPHI[ID] = cos(ARG);
    gv.SINPHI[ID] = sin(ARG);
    if (ID == gv.NPHI - 1) {
      ISIGN = 0;
    } else if (ID == gv.NPHI) {
      ISIGN = -1;
    } else if (ID == 2 * gv.NPHI - 1) {
      ISIGN = 0;
    } else if (ID == 2 * gv.NPHI) {
      ISIGN = +1;
    } else if (ID == 3 * gv.NPHI - 1) {
      ISIGN = 0;
    } else if (ID == 3 * gv.NPHI) {
      ISIGN = -1;
    }
    if (ISIGN != 0) {
      gv.S2SIGN[ID] = ISIGN;
    } else {
      gv.S2SIGN[ID] = +gv.S2SIGN[ID - 1];
    }
  }

  /* This code calculates cos and sin for discrete phi steps, omitted for brevity. */

  /* -------------------------------------------------------------------
     Scale factors
     ------------------------------------------------------------------- */
  gv.FAC = 4.0;
  if ((SIGU != 0.0) && (SIGV != 0.0))
    gv.C1 = gv.N * gv.N * FINE_STRUCTURE_CONST * BW * CUR * C_MA_A / EC / (TWOPI * SIGU * SIGV * D2 * C_M_MM * C_M_MM);

  /* -------------------------------------------------------------------
     Call the actual spectral distribution code
     ------------------------------------------------------------------- */
  spectral_distribution(IERROR, &gv);
  if (*IERROR != 0)
    return;

  /* -------------------------------------------------------------------
     Return results (E_F, SPEC0_F) to caller
     ------------------------------------------------------------------- */
  if (gv.METHOD == 4) {
    /* Dejus method => use EU[] and SPEC0[] up to NEU points */
    IP1 = 1;
    IP2 = gv.NEU;
    *NE_F = IP2 - IP1 + 1;
    for (IP = 0; IP < *NE_F; IP++) {
      E_F[IP] = gv.EU[IP];
      SPEC0_F[IP] = gv.SPEC0[IP];
    }
  } else {
    /* Walker~s methods => E array from NE1..NE2 */
    IP1 = gv.NE1;
    IP2 = gv.NE2;
    *NE_F = IP2 - IP1 + 1;
    for (IP = 0; IP < *NE_F; IP++) {
      int idx = gv.NE1 + IP;
      E_F[IP] = gv.E[idx];
      SPEC0_F[IP] = gv.SPEC0[idx];
    }
  }
  return;

  /* -------------------------------------------------------------------
     Error labels (goto targets).
     For brevity, only the main ones are shown here.
     ------------------------------------------------------------------- */
 Error900:
  printf("&USB-F-BNDERR: # of points > E_SZ=%d\n", E_SZ);
  *IERROR = -1;
  return;

 Error902:
  printf("&USB-E-INVDAT: Kx, Ky must be equal or one=0.\n");
  *IERROR = -1;
  return;

 Error905:
  printf("&USB-F-BNDERR: gv.NE > E_SZ=%d\n", E_SZ);
  *IERROR = -1;
  return;

 Error910:
  printf("&USB-E-HARMERR: No harmonics reachable in requested energy range.\n");
  *IERROR = -1;
  return;

 Error930:
  printf("&USB-E-INVDAT: Invalid method=%d\n", METHOD_F);
  *IERROR = -1;
  return;
}
