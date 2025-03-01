#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "gsl/gsl_sf_bessel.h"

#define B_SZ 1000

/* Constants */
static const double KMIN = 1.0e-3;

/**
 * @brief Computes Bessel functions of the first kind from order 0 to nmax.
 *
 * This function calculates the Bessel functions \( J_0(x) \) through \( J_{nmax}(x) \)
 * using Miller's downward recurrence algorithm. The results are stored in the
 * provided array `bs`.
 *
 * @param[in] x      The argument at which to evaluate the Bessel functions.
 * @param[in] nmax   The maximum order of the Bessel functions to compute.
 * @param[out] bs    Array to store the computed Bessel function values.
 *
 * @warning The array `bs` must be pre-allocated with at least `nmax + 1` elements.
 */
void bright_bessjn(double x, int nmax, double *bs) {
  const double BIGNO = 1.0e10;
  const double SMLNO = 1.0e-10;
  const double TWO = 2.0;
  double tox, bj, bjm, bjp, sum, sign;
  int m, i, j, jsum;

  tox = TWO / fabs(x);
  /* Make m even */
  m = 2 * ((nmax + 1) / 2);
  jsum = 0;
  sum = 0.0;
  bjp = 0.0;
  bj = 1.0;

  /* Miller's downward recurrence */
  for (j = m; j >= 1; j--) {
    bjm = j * tox * bj - bjp;
    bjp = bj;
    bj = bjm;
    if (fabs(bj) > BIGNO) {
      bj *= SMLNO;
      bjp *= SMLNO;
      sum *= SMLNO;
      for (i = j; i <= nmax; i++) {
        bs[i] *= SMLNO;
      }
    }
    if (jsum != 0)
      sum += bj;
    jsum = 1 - jsum;
    if (j <= nmax)
      bs[j] = bjp;
  }
  bs[0] = bj;
  sum = 2.0 * sum - bj;

  /* Normalize */
  if (x < 0.0) {
    sign = -1.0;
    for (j = 0; j <= nmax; j++) {
      sign = -sign;
      bs[j] = sign * bs[j] / sum;
    }
  } else {
    for (j = 0; j <= nmax; j++) {
      bs[j] /= sum;
    }
  }
}

/**
 * @brief Computes Bessel functions of the first kind using the GNU Scientific Library (GSL).
 *
 * This function serves as a replacement for `bright_bessjn` by utilizing GSL's
 * `gsl_sf_bessel_Jn_array` function to compute \( J_0(x) \) through \( J_{nmax}(x) \).
 * It handles negative values of `x` by applying the identity \( J_n(-x) = (-1)^n J_n(x) \).
 *
 * @param[in] x      The argument at which to evaluate the Bessel functions.
 * @param[in] nmax   The maximum order of the Bessel functions to compute.
 * @param[out] bs    Array to store the computed Bessel function values.
 *
 * @note Requires the GSL library. Ensure that `bs` is pre-allocated with at least `nmax + 1` elements.
 */
void bright_bessjn_gsl(double x, int nmax, double *bs) {
  /* GSL provides J0, J1, ..., Jnmax in one call. */
  /* If x < 0, use Jn(-x)=(-1)^n Jn(x). */
  int n;
  if (nmax < 0) return; /* Safety check */

  if (x >= 0.0) {
    /* Directly compute J0..Jnmax for x>=0 */
    gsl_sf_bessel_Jn_array(0, nmax, x, bs);
  } else {
    /* For negative x, compute for |x|, then flip sign for odd n. */
    double pos_x = -x;
    double *temp = (double *)malloc((nmax + 1) * sizeof(double));
    if (!temp) {
      fprintf(stderr, "Memory allocation failed\n");
      exit(EXIT_FAILURE);
    }
    gsl_sf_bessel_Jn_array(0, nmax, pos_x, temp);
    for (n = 0; n <= nmax; n++) {
      /* Jn(-x) = (-1)^n * Jn(x) */
      if (n % 2 == 0)
        bs[n] = temp[n];  /* even n => + */
      else
        bs[n] = -temp[n]; /* odd n => - */
    }
    free(temp);
  }
}

/**
 * @brief Computes Stokes parameters for a regular or flipped plane device.
 *
 * This subroutine, translated from Fortran to C, calculates the Stokes parameters
 * \( s0, s1, s2, s3 \) based on the input indices and physical parameters. It utilizes
 * Bessel functions computed via `bright_bessjn` to perform the calculations.
 *
 * @param[in]  i        The harmonic index.
 * @param[in]  k_val    The wave number value.
 * @param[in]  alpha    The alpha parameter.
 * @param[in]  cosphi   The cosine of the phi angle.
 * @param[in]  sinphi   The sine of the phi angle.
 * @param[out] s0       Pointer to store the computed Stokes parameter S0.
 * @param[out] s1       Pointer to store the computed Stokes parameter S1.
 * @param[out] s2       Pointer to store the computed Stokes parameter S2.
 * @param[out] s3       Pointer to store the computed Stokes parameter S3.
 *
 */
void bright1(int i, double k_val, double alpha, double cosphi, double sinphi,
             double *s0, double *s1, double *s2, double *s3) {
  /* Array for Bessel values */
  double axr, ayr;
  double bx[B_SZ + 1], by[B_SZ + 1];
  int n, nx, ny, np, nm, im;
  double alpha2, k2, a, c, x, y, absx;
  double sum1, sum2, sum3, sign;
  const double EPS = 1e-6;
  const double EPSX = 1e-5;
  const double HALF = 0.5;
  const double TWO = 2.0;
  int lodd = (i % 2 != 0); /* .true. if i is odd */

  alpha2 = alpha * alpha;
  k2 = k_val * k_val;
  a = 1.0 + 0.5 * k2 + alpha2;
  c = -k_val * alpha * cosphi;
  x = (TWO * i / a) * (-c);
  y = (0.25 * i / a) * (k2);
  absx = fabs(x);

  /* Find Bessel function ranges */
  nx = (int)(6.20 + 1.41 * absx + 1.0);
  if (nx > B_SZ) {
    nx = B_SZ;
    printf(" &bright1-W-DIMERR, Dimension error\n");
    printf(" Array of Bessel function values out of bound for x\n");
    printf(" - upper bound reset to %d\n", B_SZ);
  }
  ny = (int)(6.20 + 1.41 * fabs(y) + 1.0);
  if (ny > B_SZ) {
    ny = B_SZ;
    printf(" &bright1-W-DIMERR, Dimension error\n");
    printf(" Array of Bessel function values out of bound for y\n");
    printf(" - upper bound reset to %d\n", B_SZ);
  }

  /* Clear arrays */
  for (n = 0; n <= B_SZ; n++) {
    bx[n] = 0.0;
    by[n] = 0.0;
  }

  if (absx >= EPSX) {
    bright_bessjn(x, nx, bx);
  }
  bright_bessjn(fabs(y), ny, by);

  /* Adjust ny if by[ny] < EPS */
  while (ny > 0 && by[ny] < EPS) {
    ny--;
  }

  sum1 = 0.0;
  sum2 = 0.0;
  sum3 = 0.0;

  if (absx < EPSX) {
    /* For very small x, handle odd/even harmonic separately */
    if (lodd) {
      /* Odd harmonic (e.g. i=1,3,5,...) */
      sum1 = 0.0;
      sum3 = 0.0;
      np = -(-i + 1) / 2;
      nm = -(-i - 1) / 2;

      if (np <= ny && np >= 0) {
        if (np % 2 != 0) {
          sum1 = -by[np];
          sum3 = -by[np];
        } else {
          sum1 = by[np];
          sum3 = by[np];
        }
      }
      if (nm <= ny && nm >= 0) {
        if (nm % 2 != 0) {
          sum1 += by[nm];
          sum3 -= by[nm];
        } else {
          sum1 -= by[nm];
          sum3 += by[nm];
        }
      }
      sum1 = HALF * x * sum1;
    } else {
      /* Even harmonic (e.g. i=2,4,6,...) */
      sum1 = 0.0;
      nm = i / 2;
      if (nm <= ny && nm >= 0) {
        if (nm % 2 != 0)
          sum1 = -by[nm];
        else
          sum1 = by[nm];
      }

      sum3 = 0.0;
      np = -(-i + 2) / 2;
      nm = -(-i - 2) / 2;
      if (np <= ny && np >= 0) {
        if (np % 2 != 0)
          sum3 = -by[np];
        else
          sum3 = by[np];
      }
      if (nm <= ny && nm >= 0) {
        if (nm % 2 != 0)
          sum3 += by[nm];
        else
          sum3 -= by[nm];
      }
      sum3 = HALF * x * sum3;
    }
  } else {
    /* Larger x, do full expansions */
    sign = 1.0;
    if (i <= nx)
      sum1 = by[0] * bx[i];
    for (n = 1; n <= ny; n++) {
      sign = -sign;
      int np_ = 2 * n + i;
      int nm_ = i - 2 * n;
      im = abs(nm_);

      double bxp = (np_ <= nx && np_ >= 0) ? bx[np_] : 0.0;
      double bxm = 0.0;
      if (im <= nx && im >= 0) {
        if (lodd && (nm_ < 0))
          bxm = -bx[im];
        else
          bxm = bx[im];
      }

      double bx1 = bxp + sign * bxm;
      double bx2 = bxp - sign * bxm;
      sum1 += by[n] * bx1;
      sum2 += n * by[n] * bx2;
    }
    sum3 = TWO / x * (i * sum1 + 2.0 * sum2);
  }

  axr = (2.0 * i / a) * (alpha * cosphi * sum1 - 0.5 * k_val * sum3);
  ayr = (2.0 * i / a) * (alpha * sinphi * sum1);

  /* Stokes */
  *s0 = axr * axr + ayr * ayr;
  *s1 = axr * axr - ayr * ayr;
  *s2 = 2.0 * axr * ayr;
  *s3 = 0.0;
}

/**
 * @brief Computes Stokes parameters for an elliptical device.
 *
 * This subroutine, translated from Fortran to C, calculates the Stokes parameters
 * \( s0, s1, s2, s3 \) for an elliptical device based on the input indices and
 * physical parameters. It utilizes Bessel functions computed via `bright_bessjn` to
 * perform the calculations and handles both small and large argument scenarios.
 *
 * @param[in]  i        The harmonic index.
 * @param[in]  kx_val   The wave number value in the x-direction.
 * @param[in]  ky_val   The wave number value in the y-direction.
 * @param[in]  alpha    The alpha parameter.
 * @param[in]  cosphi   The cosine of the phi angle.
 * @param[in]  sinphi   The sine of the phi angle.
 * @param[out] s0       Pointer to store the computed Stokes parameter S0.
 * @param[out] s1       Pointer to store the computed Stokes parameter S1.
 * @param[out] s2       Pointer to store the computed Stokes parameter S2.
 * @param[out] s3       Pointer to store the computed Stokes parameter S3.
 *
 * @note This function handles various cases based on the magnitude of `x` and `y`.
 *       It includes complex calculations involving trigonometric recurrences and
 *       Bessel function combinations.
 */
void bright3(int i, double kx_val, double ky_val, double alpha,
             double cosphi, double sinphi,
             double *s0, double *s1, double *s2, double *s3) {
  double axr, ayr, axi, ayi;
  double bx[B_SZ + 1], by[B_SZ + 1];
  double alpha2, k2, a, x, y, absy, phi;
  double co, so, cp0, sp0, cm0, sm0;
  double cp1p, sp1p, cm1p, sm1p;
  double cp1m, sp1m, cm1m, sm1m;
  double sum0r, sum0i, sum1pr, sum1pi, sum1mr, sum1mi;
  double sign;
  int nx, ny, n;
  int np0, np1p, np1m, nm0, nm1p, nm1m, im0, im1p, im1m;
  const double EPS = 1e-6;
  const double EPSX = 1e-5;
  const double EPSY = 1e-5;
  const double TWO = 2.0;
  int lodd = (i % 2 != 0); /* .true. if i is odd */

  alpha2 = alpha * alpha;
  k2 = kx_val * kx_val + ky_val * ky_val;
  a = 1.0 + 0.5 * k2 + alpha2;

  /* x, y definitions from Fortran code */
  x = (TWO * i / a) * alpha * sqrt(kx_val * kx_val * sinphi * sinphi + ky_val * ky_val * cosphi * cosphi);
  y = (0.25 * i / a) * (ky_val * ky_val - kx_val * kx_val);
  absy = fabs(y);
  phi = atan2(kx_val * sinphi, ky_val * cosphi);

  co = cos(2.0 * phi);
  so = sin(2.0 * phi);
  cp0 = cos(i * phi);
  sp0 = sin(i * phi);
  cm0 = cp0; /* in original code cm0=cp0, sm0=sp0 */
  sm0 = sp0;
  cp1p = cos((i + 1) * phi);
  sp1p = sin((i + 1) * phi);
  cm1p = cp1p;
  sm1p = sp1p;
  cp1m = cos((i - 1) * phi);
  sp1m = sin((i - 1) * phi);
  cm1m = cp1m;
  sm1m = sp1m;

  sum0r = 0.0;
  sum0i = 0.0;
  sum1pr = 0.0;
  sum1pi = 0.0;
  sum1mr = 0.0;
  sum1mi = 0.0;

  /* Determine Bessel ranges */
  nx = (int)(6.20 + 1.41 * x + 1.0);
  if (nx > B_SZ) {
    nx = B_SZ;
    printf(" &bright3-W-DIMERR, Dimension error\n");
    printf(" Array of Bessel function values out of bound for x\n");
    printf(" - upper bound reset to %d\n", B_SZ);
  }
  ny = (int)(6.20 + 1.41 * absy + 1.0);
  if (ny > B_SZ) {
    ny = B_SZ;
    printf(" &bright3-W-DIMERR, Dimension error\n");
    printf(" Array of Bessel function values out of bound for y\n");
    printf(" - upper bound reset to %d\n", B_SZ);
  }

  /* Clear arrays */
  for (n = 0; n <= B_SZ; n++) {
    bx[n] = 0.0;
    by[n] = 0.0;
  }

  if (x >= EPSX) {
    bright_bessjn(x, nx, bx);
  }
  if (absy >= EPSY) {
    bright_bessjn(absy, ny, by);
    while (ny > 0 && by[ny] < EPS) {
      ny--;
    }
  }

  /* Various cases for small/large x, y */
  if (x < EPSX && absy < EPSY) {
    /* Both x, y small */
    if (i == 1) {
      sum0r = cos(phi) * 0.5 * x;
      sum0i = -sin(phi) * 0.5 * x;
      sum1mr = 1.0;
    } else if (i == 2) {
      sum1mr = cos(phi) * 0.5 * x;
      sum1mi = -sin(phi) * 0.5 * x;
    }
  } else if (x < EPSX) {
    /* x small */
    if (lodd) {
      int np = -(-i + 1) / 2;
      int nm = -(-i - 1) / 2;
      if (np <= ny && np >= 0) {
        if (np % 2 != 0) {
          sum0r = -by[np];
          sum0i = -by[np];
          sum1mr = -by[np];
        } else {
          sum0r = by[np];
          sum0i = by[np];
          sum1mr = by[np];
        }
      }
      if (nm <= ny && nm >= 0) {
        if (nm % 2 != 0) {
          sum0r += by[nm];
          sum0i -= by[nm];
          sum1pr = -by[nm];
        } else {
          sum0r -= by[nm];
          sum0i += by[nm];
          sum1pr = by[nm];
        }
      }
      sum0r = cos(phi) * 0.5 * x * sum0r;
      sum0i = -sin(phi) * 0.5 * x * sum0i;
    } else {
      int nm_ = i / 2;
      double sum0temp = 0.0;
      if (nm_ <= ny && nm_ >= 0) {
        if (nm_ % 2 != 0)
          sum0temp = -by[nm_];
        else
          sum0temp = by[nm_];
      }
      int np = -(-i + 2) / 2;
      nm_ = -(-i - 2) / 2;
      double sum1p = 0.0, sum1m = 0.0;
      if (np <= ny && np >= 0) {
        if (np % 2 != 0)
          sum1m = -by[np];
        else
          sum1m = by[np];
      }
      if (nm_ <= ny && nm_ >= 0) {
        if (nm_ % 2 != 0)
          sum1p = -by[nm_];
        else
          sum1p = by[nm_];
      }
      sum1pr = cos(phi) * 0.5 * x * (sum0temp - sum1p);
      sum1pi = -sin(phi) * 0.5 * x * (sum0temp + sum1p);
      sum1mr = cos(phi) * 0.5 * x * (sum1m - sum0temp);
      sum1mi = -sin(phi) * 0.5 * x * (sum1m + sum0temp);
    }
  } else if (absy < EPSY) {
    /* y small */
    if (i <= nx) {
      sum0r = cos(i * phi) * bx[i];
      sum0i = -sin(i * phi) * bx[i];
    }
    if ((i + 1) <= nx) {
      sum1pr = cos((i + 1) * phi) * bx[i + 1];
      sum1pi = -sin((i + 1) * phi) * bx[i + 1];
    }
    if ((i - 1) >= 0 && (i - 1) <= nx) {
      sum1mr = cos((i - 1) * phi) * bx[i - 1];
      sum1mi = -sin((i - 1) * phi) * bx[i - 1];
    }
  } else {
    /* Large x,y: full expansions */
    if (i <= nx) {
      sum0r = cos(i * phi) * by[0] * bx[i];
      sum0i = -sin(i * phi) * by[0] * bx[i];
    }
    if ((i + 1) <= nx) {
      sum1pr = cos((i + 1) * phi) * by[0] * bx[i + 1];
      sum1pi = -sin((i + 1) * phi) * by[0] * bx[i + 1];
    }
    if ((i - 1) >= 0 && (i - 1) <= nx) {
      sum1mr = cos((i - 1) * phi) * by[0] * bx[i - 1];
      sum1mi = -sin((i - 1) * phi) * by[0] * bx[i - 1];
    }
    sign = 1.0;
    for (n = 1; n <= ny; n++) {
      sign = -sign;
      np0 = 2 * n + i;
      np1p = 2 * n + i + 1;
      np1m = 2 * n + i - 1;
      nm0 = i - 2 * n;
      nm1p = i - 2 * n + 1;
      nm1m = i - 2 * n - 1;
      im0 = abs(nm0);
      im1p = abs(nm1p);
      im1m = abs(nm1m);

      /* --- Replicate the trig recurrence from Fortran exactly --- */
      {
        double cps, cms;

        /* Update (cp0, sp0) via co, so */
        cps = cp0;
        cms = cm0;
        cp0 = cp0 * co - sp0 * so; /* cp0 = cp0 * co - sp0 * so */
        cm0 = cm0 * co + sm0 * so;
        sp0 = sp0 * co + cps * so;
        sm0 = sm0 * co - cms * so;

        /* Update (cp1p, sp1p) */
        cps = cp1p;
        cms = cm1p;
        cp1p = cp1p * co - sp1p * so;
        cm1p = cm1p * co + sm1p * so;
        sp1p = sp1p * co + cps * so;
        sm1p = sm1p * co - cms * so;

        /* Update (cp1m, sp1m) */
        cps = cp1m;
        cms = cm1m;
        cp1m = cp1m * co - sp1m * so;
        cm1m = cm1m * co + sm1m * so;
        sp1m = sp1m * co + cps * so;
        sm1m = sm1m * co - cms * so;
      }
      /* Now get bx(...) for positive/ negative order, as Fortran does. */
      /* bxp0r = cp0*bx(np0), bxp0i= sp0*bx(np0), etc. if np0 <= nx */
      double bxp0r = 0.0, bxp0i = 0.0;
      double bxp1pr = 0.0, bxp1pi = 0.0;
      double bxp1mr = 0.0, bxp1mi = 0.0;
      double bxm0r = 0.0, bxm0i = 0.0;
      double bxm1pr = 0.0, bxm1pi = 0.0;
      double bxm1mr = 0.0, bxm1mi = 0.0;

      if (np0 >= 0 && np0 <= nx) {
        bxp0r = cp0 * bx[np0]; /* real part = cos(...) * J_{np0} */
        bxp0i = sp0 * bx[np0]; /* imag part = sin(...) * J_{np0} */
      }
      if (np1p >= 0 && np1p <= nx) {
        bxp1pr = cp1p * bx[np1p];
        bxp1pi = sp1p * bx[np1p];
      }
      if (np1m >= 0 && np1m <= nx) {
        bxp1mr = cp1m * bx[np1m];
        bxp1mi = sp1m * bx[np1m];
      }

      /* For negative order (nm0, etc.) in Fortran, sign flips if odd/even, etc. */
      if (im0 <= nx) {
        double factor = 1.0;
        /* if harmonic is odd and nm0<0 => factor=-1 (Fortran code) */
        if (lodd && (nm0 < 0)) {
          factor = -1.0;
        }
        bxm0r = factor * cm0 * bx[im0];
        bxm0i = factor * sm0 * bx[im0];
      }
      if (im1p <= nx) {
        double factor = 1.0;
        /* if harmonic is even and nm1p<0 => factor=-1 in Fortran code */
        if ((!lodd) && (nm1p < 0)) {
          factor = -1.0;
        }
        bxm1pr = factor * cm1p * bx[im1p];
        bxm1pi = factor * sm1p * bx[im1p];
      }
      if (im1m <= nx) {
        double factor = 1.0;
        if ((!lodd) && (nm1m < 0)) {
          factor = -1.0;
        }
        bxm1mr = factor * cm1m * bx[im1m];
        bxm1mi = factor * sm1m * bx[im1m];
      }

      /* Combine per the Fortran code: bxp + sign*bxm. */
      {
        double bx10r = bxp0r + sign * bxm0r;
        double bx10i = bxp0i + sign * bxm0i;
        double bx11pr = bxp1pr + sign * bxm1pr;
        double bx11pi = bxp1pi + sign * bxm1pi;
        double bx11mr = bxp1mr + sign * bxm1mr;
        double bx11mi = bxp1mi + sign * bxm1mi;

        /* Then accumulate sums (noting Fortran uses sum0r -= by[n]*bx10i for imaginary) */
        sum0r += by[n] * bx10r;
        sum0i -= by[n] * bx10i;
        sum1pr += by[n] * bx11pr;
        sum1pi -= by[n] * bx11pi;
        sum1mr += by[n] * bx11mr;
        sum1mi -= by[n] * bx11mi;
      }
    } /* end for n=1..ny */
  }

  {
    double denom = (double)i / a;
    double sum1pr_mr_r = sum1pr + sum1mr; /* real parts for horizontal pol. */
    double sum1pi_mi_r = sum1pi + sum1mi; /* imaginary parts for horizontal */
    double sum1pi_mi_d = sum1pi - sum1mi; /* difference for vertical pol.   */
    double sum1pr_mr_d = sum1pr - sum1mr; /* difference for vertical pol.   */

    axr = denom * (2.0 * alpha * cosphi * sum0r - ky_val * sum1pr_mr_r);
    axi = denom * (2.0 * alpha * cosphi * sum0i - ky_val * sum1pi_mi_r);
    ayr = denom * (2.0 * alpha * sinphi * sum0r + kx_val * sum1pi_mi_d);
    ayi = denom * (2.0 * alpha * sinphi * sum0i - kx_val * sum1pr_mr_d);
  }
  {
    double axr2 = axr * axr;
    double axi2 = axi * axi;
    double ayr2 = ayr * ayr;
    double ayi2 = ayi * ayi;

    *s0 = axr2 + axi2 + ayr2 + ayi2;
    *s1 = axr2 + axi2 - ayr2 - ayi2;
    *s2 = 2.0 * (axr * ayr + axi * ayi);
    *s3 = 2.0 * (axi * ayr - axr * ayi);
  }
}

/* ----------------------------------------------------------------------
   brighte subroutine (Fortran -> C)
   ---------------------------------------------------------------------- */
void brighte(int i, double alpha, double cosphi, double sinphi,
             double *s0, double *s1, double *s2, double *s3, double kx, double ky) {
  /* Checks correspond to the Fortran if statements */
  if ((kx < KMIN) && (ky > KMIN)) {
    /* Regular plane device */
    bright1(i, ky, alpha, cosphi, sinphi, s0, s1, s2, s3);
  } else if ((kx > KMIN) && (ky < KMIN)) {
    /* Flipped plane device */
    bright1(i, kx, alpha, sinphi, -cosphi, s0, s1, s2, s3);
    *s1 = -(*s1);
    *s2 = -(*s2);
  } else {
    /* Elliptical device */
    bright3(i, kx, ky, alpha, cosphi, sinphi,
            s0, s1, s2, s3);
  }
}
