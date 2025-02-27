/************************************************************************* \
* Copyright (c) 2019 The University of Chicago, as Operator of Argonne
* National Laboratory.
* This file is distributed subject to a Software License Agreement found
* in the file LICENSE that is included with this distribution. 
\*************************************************************************/

/* file: lorentzian.c
 * purpose: compute electric field due to lorentzian charge distribution
 *
 * The charge distribution is Q*a*b/[\pi^2 (a^2 + x^2) (b^2 + y^2)]
 *
 * Ryan Lindberg (2019).
 */

#include "mdb.h"
#include "track.h"

#define TWOPI PIx2

/* Gives V/Q, the field*length (Volts) per Coulomb */

void evaluateVoltageFromLorentzian(double *Eperp, double a, double b, double x, double y) {
  /* what I found to be convenient auxilliary parameters */
  double bMINx, bPLUSx, aMINy, aPLUSy;
  double denomBM_AM, denomBP_AM, denomBM_AP, denomBP_AP;
  double denom, tempx, tempy, aCubed, bCubed;
  double epsx, epsy;
  double arcTANa_b, arcTANb_a;

  /* Check to see if general expression can be used; it has numerical problems when |x| -> |b| and |y| -> |a| */
  if ((fabs(fabs(x / b) - 1) < 0.001) && (fabs(fabs(y / a) - 1) < 0.001)) {
    /* Use expansion accurate to better than 3x10^{-10} and make calculation assuming x > 0, y > 0 */
    arcTANa_b = atan2(a, b);
    arcTANb_a = atan2(b, a);

    epsx = fabs(x) - b;
    epsy = fabs(y) - a;
    denom = 1.0 / (a * a + b * b);
    tempx = (-2.0 * a * b + PI * (a * a + 2.0 * b * b)) * denom;
    Eperp[0] = (0.5 * tempx - arcTANa_b) / b;
    Eperp[1] = (0.5 * tempx - arcTANb_a) / a;

    aCubed = a * a * a;
    bCubed = b * b * b;
    tempx = 2.0 * arcTANa_b / denom;
    tempy = 2.0 * arcTANb_a / denom;
    denom = denom * denom;
    tempx += a * (2.0 * (aCubed * a * b + bCubed * (4.0 * a * a - b * b)) - PI * (aCubed * (a * a + b * b) + 2.0 * a * b * bCubed)) * denom;
    tempy += b * (2.0 * (bCubed * b * a + aCubed * (4.0 * b * b - a * a)) - PI * (bCubed * (b * b + a * a) + 2.0 * b * a * aCubed)) * denom;

    Eperp[0] += 0.25 * epsx * tempx / (a * a * b * b) - epsy * 0.5 * PI * a * b * denom;
    Eperp[1] += 0.25 * epsy * tempy / (a * a * b * b) - epsx * 0.5 * PI * a * b * denom;

    denom = denom / (a * a + b * b);
    tempx = (3.0 * a * b * b - aCubed) * (6.0 * a * a * b + 2.0 * bCubed + 3.0 * PI * aCubed);
    tempx = (tempx * denom - 6.0 * arcTANa_b) / (12.0 * aCubed);
    tempy = (3.0 * b * a * a - bCubed) * (6.0 * b * b * a + 2.0 * aCubed + 3.0 * PI * bCubed);
    tempy = (tempy * denom - 6.0 * arcTANb_a) / (12.0 * bCubed);

    Eperp[0] += epsx * epsy * tempx;
    Eperp[1] += epsx * epsy * tempy;

    tempx = -aCubed * (6.0 * a * a * b + 16.0 * bCubed) + 3.0 * PI * (aCubed * (aCubed + 3.0 * a * b * b) + 2.0 * bCubed * bCubed) - 90.0 * a * b * b * bCubed;
    tempx = (tempx * denom - 6.0 * arcTANa_b) / (24.0 * bCubed);
    tempy = -bCubed * (6.0 * b * b * a + 16.0 * aCubed) + 3.0 * PI * (bCubed * (bCubed + 3.0 * b * a * a) + 2.0 * aCubed * aCubed) - 90.0 * b * a * a * aCubed;
    tempy = (tempy * denom - 6.0 * arcTANb_a) / (24.0 * aCubed);

    Eperp[0] += epsx * epsx * tempx - epsy * epsy * (tempx + 4.0 * a * b * b * denom);
    Eperp[1] += epsy * epsy * tempy - epsx * epsx * (tempy + 4.0 * b * a * a * denom);
    /* reflect E-field if x or y is negative */
    if (x < 0)
      Eperp[0] *= -1.0;
    if (y < 0)
      Eperp[1] *= -1.0;
  } else { /* Using the full expression */
    bMINx = b - x;
    bPLUSx = b + x;
    aMINy = a - y;
    aPLUSy = a + y;
    denomBM_AM = 1.0 / (bMINx * bMINx + aMINy * aMINy);
    denomBP_AM = 1.0 / (bPLUSx * bPLUSx + aMINy * aMINy);
    denomBM_AP = 1.0 / (bMINx * bMINx + aPLUSy * aPLUSy);
    denomBP_AP = 1.0 / (bPLUSx * bPLUSx + aPLUSy * aPLUSy);

    Eperp[0] = (bMINx * (denomBM_AM + denomBM_AP) + bPLUSx * (denomBP_AM + denomBP_AP)) * atan2(x, a);
    Eperp[0] += (-bMINx * (denomBM_AM + denomBM_AP) + bPLUSx * (denomBP_AM + denomBP_AP)) * 0.5 * PI;

    Eperp[1] = (aMINy * (denomBM_AM + denomBP_AM) + aPLUSy * (denomBM_AP + denomBP_AP)) * atan2(y, b);
    Eperp[1] += (-aMINy * (denomBM_AM + denomBP_AM) + aPLUSy * (denomBM_AP + denomBP_AP)) * 0.5 * PI;

    denom = denomBM_AM * denomBP_AM * denomBM_AP * denomBP_AP;
    tempx = 8.0 * a * y * x * (bMINx * bPLUSx * (2.0 * (a * a + y * y) + 3.0 * b * b + x * x) - aMINy * aMINy * aPLUSy * aPLUSy) * atan2(y, b);
    tempx -= 4.0 * a * b * x * (aMINy * aPLUSy * (2.0 * (b * b + x * x) + 3.0 * y * y + a * a) + bPLUSx * bPLUSx * bMINx * bMINx) * log((a * a + x * x) / (b * b + y * y));

    tempy = 8.0 * b * x * y * (aMINy * aPLUSy * (2.0 * (b * b + x * x) + 3.0 * a * a + y * y) - bMINx * bMINx * bPLUSx * bPLUSx) * atan2(x, a);
    tempy += 4.0 * b * a * y * (bMINx * bPLUSx * (2.0 * (a * a + y * y) + 3.0 * x * x + b * b) + aPLUSy * aPLUSy * aMINy * aMINy) * log((a * a + x * x) / (b * b + y * y));

    Eperp[0] += tempx * denom;
    Eperp[1] += tempy * denom;
  }

  Eperp[0] *= 1. / (4 * sqr(PI) * epsilon_o);
  Eperp[1] *= 1. / (4 * sqr(PI) * epsilon_o);

  return;
}
