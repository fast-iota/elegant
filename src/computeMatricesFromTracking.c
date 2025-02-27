/* Copyright 1994, 2015 by Michael Borland and Argonne National Laboratory,
 * all rights reserved.
 */
/* routine: computeMatricesFromTracking()
 * purpose: compute R, T, and Q Transport matrices from tracking data.
 *
 * This version has been modified to correctly subtract off certain high-order
 * terms (i.e., >cubic) after fitting, rather than just subtracting off
 * terms included in the transport matrix expansion.  The terms subtracted off
 * are those that depend only on one of the initial coordinates.  While this 
 * won't change the matrices, it will change the residuals. 
 *
 * Michael Borland, 1988, 1989, 1990, 2015
 */
#include "mdb.h"
#include "track.h"

#define CONSTRAIN_EQ 0
#define CONSTRAIN_NEQ 1
#define UNCONSTRAINED 2

#define DEBUG1 0
#define DEBUG2 0
#define DEBUG3 0

#ifndef HUGE_VAL
#  define HUGE_VAL HUGE
#endif

typedef struct {
  double *coefficient, *coefficient_error;
  int order_of_fit; /* 1=linear, 2=quadratic...*/
  double chi_squared;
} FIT;

int selectInitialVectors(double **set_i, double **set_f, double **set_e,
                         int n_max_select, double **initial, double **final, double **error,
                         int n_vectors, double *resolution, double *constraint, int *constraint_code);

void copyFit(FIT *target, FIT *source);
void showSelectedPairs(FILE *fp, double **ini, double **fin, double **err, int n);
void checkAssignedErrors(FILE *fp, double **ini, double **fin, double **err, int n);
int findBestFit(FIT *best_fit, double **set_i, double **set_f, double **set_error,
                int n_pairs, int i, int j, int type, int max_order);

VMATRIX *computeMatricesFromTracking(
  FILE *fpo_ma,
  double **initial, double **final, double **error,
  double *step_size, double *maximum_value,
  int n_points1, int n_points_total,
  int max_order_of_fits,
  int verbose) {
  long i, j, k, l, m;
  long n_points_max, n_pairs1, n_pairs2;
  double sum, term;
  double **set1_i, **set1_f, **set1_error;
  double **set2_i, **set2_f, **set2_error;
  double *resolution, *constraint, Tijk;
  int *constraint_code;
  FIT best_fit, best_fit1, best_fit2;
  FIT **saved_fit;
  long *n_fits;
  VMATRIX *M;
  double *C, **R, ***T, ****Q;

  M = tmalloc(sizeof(*M));
  M->order = 3;
  initialize_matrices(M, M->order);
  null_matrices(M, 0);
  C = M->C;
  R = M->R;
  T = M->T;
  Q = M->Q;

  best_fit.coefficient = best_fit.coefficient_error = NULL;
  best_fit1.coefficient = best_fit1.coefficient_error = NULL;
  best_fit2.coefficient = best_fit2.coefficient_error = NULL;

  n_points_max = 0;
  resolution = tmalloc(sizeof(*resolution) * 6);
  constraint = tmalloc(sizeof(*constraint) * 6);
  n_fits = tmalloc(sizeof(*n_fits) * 6);
  constraint_code = tmalloc(sizeof(*constraint_code) * 6);
  n_points_max = n_points1;
  for (i = 5; i >= 0; i--) {
    n_fits[i] = C[i] = 0;
    resolution[i] = step_size[i] / 10.0;
    R[i][i] = 1;
  }
  /* allocate parallel arrays for storing sets of initial and final
   * coordinate pairs. 
   */
  set1_i = (double **)zarray_2d(sizeof(**set1_i), n_points_max, 6);
  set1_f = (double **)zarray_2d(sizeof(**set1_f), n_points_max, 6);
  set1_error = (double **)zarray_2d(sizeof(**set1_error), n_points_max, 6);
  set2_i = (double **)zarray_2d(sizeof(**set2_i), n_points_max, 6);
  set2_f = (double **)zarray_2d(sizeof(**set2_f), n_points_max, 6);
  set2_error = (double **)zarray_2d(sizeof(**set2_error), n_points_max, 6);

  /* Compute Ci, Rij, Tijj, and Qijjj. */
  saved_fit = (FIT **)zarray_2d(sizeof(**saved_fit), 6, 6);
  for (j = 5; j >= 0; j--) {
    /* select a set of initial/final vector pairs such that initial 
     * coordinates are zero for coodinates other than the jth.
     */
    for (i = 5; i >= 0; i--)
      if (i != j) {
        constraint[i] = 0.0;
        constraint_code[i] = CONSTRAIN_EQ;
      } else {
        constraint_code[i] = UNCONSTRAINED;
      }

    if ((n_pairs1 = selectInitialVectors(set1_i, set1_f, set1_error,
                                         n_points_max, initial, final, error,
                                         n_points_total, resolution, constraint,
                                         constraint_code)) < 3) {
      fprintf(fpo_ma, "too few vector pairs for computing R[i][%ld]\n",
              j + 1);
      for (i = 5; i >= 0; i--)
        saved_fit[i][j].order_of_fit = 0;
      continue;
    }
    if (verbose > 1)
      fprintf(fpo_ma, "%ld pairs selected for R[i][%ld]\n",
              n_pairs1, j + 1);
#if DEBUG1
    if (verbose <= 1)
      fprintf(fpo_ma, "%ld pairs selected for R[i][%ld]\n",
              n_pairs1, j + 1);
    showSelectedPairs(fpo_ma, set1_i, set1_f, set1_error, n_pairs1);
#endif

    /* use the vector pairs to get the matrix elements Rij, 
     * Tijj, and Qijjj for i=1 to 6.  Also save the constant terms from
     * the fits.
     */
    for (i = 5; i >= 0; i--) {
      if (!findBestFit(&best_fit, set1_i, set1_f, set1_error,
                       n_pairs1, i, j, 10 * (i * 10 + j) + 1, max_order_of_fits)) {
        fprintf(fpo_ma, "unable to find fit for C[%ld], R[%ld][%ld], T[%ld][%ld][%ld], and Q[%ld][%ld][%ld][%ld]\n",
                i + 1, i + 1, j + 1, i + 1, j + 1, j + 1, i + 1, j + 1, j + 1, j + 1);
        continue;
      }
      copyFit(saved_fit[i] + j, &best_fit);
      C[i] += best_fit.coefficient[0];
      n_fits[i] += 1;
      R[i][j] = best_fit.coefficient[1];
      if (verbose > 1) {
        fprintf(fpo_ma, "order of fit for R[%ld][%ld], T[%ld][%ld][%ld], Q[%ld][%ld][%ld][%ld] = %d\n",
                i + 1, j + 1,
                i + 1, j + 1, j + 1,
                i + 1, j + 1, j + 1, j + 1,
                best_fit.order_of_fit);
        fprintf(fpo_ma, "    R[%ld][%ld]      = %.16le +/- %.16le\n",
                i + 1, j + 1, R[i][j],
                best_fit.coefficient_error[1]);
      }
#if DEBUG1
      if (verbose <= 1) {
        fprintf(fpo_ma, "order of fit for R[%ld][%ld], T[%ld][%ld][%ld], Q[%ld][%ld][%ld][%ld] = %d\n",
                i + 1, j + 1,
                i + 1, j + 1, j + 1,
                i + 1, j + 1, j + 1, j + 1,
                best_fit.order_of_fit);
        fprintf(fpo_ma, "    R[%ld][%ld]      = %.16le +/- %.16le\n",
                i + 1, j + 1, R[i][j],
                best_fit.coefficient_error[1]);
      }
      fprintf(fpo_ma, "chi-squared = %le \n", best_fit.chi_squared);
      fprintf(fpo_ma, "coefficients of fit:\n");
      for (k = 0; k < best_fit.order_of_fit + 1; k++)
        fprintf(fpo_ma, "a[%ld] = %le +/- %le\n",
                k, best_fit.coefficient[k], best_fit.coefficient_error[k]);
      fflush(fpo_ma);
#endif
      if (best_fit.order_of_fit >= 2) {
        T[i][j][j] = best_fit.coefficient[2];
        if (verbose > 1)
          fprintf(fpo_ma, "    T[%ld][%ld][%ld] = %.16le +/- %.16le\n",
                  i + 1, j + 1, j + 1, T[i][j][j],
                  best_fit.coefficient_error[2]);
      }
      if (best_fit.order_of_fit >= 3) {
        Q[i][j][j][j] = best_fit.coefficient[3];
        if (verbose > 1)
          fprintf(fpo_ma, "    Q[%ld][%ld][%ld][%ld] = %.16le +/- %.16le\n",
                  i + 1, j + 1, j + 1, j + 1, Q[i][j][j][j],
                  best_fit.coefficient_error[3]);
      }
      if (best_fit.order_of_fit >= 4 && verbose > 1) {
        fprintf(fpo_ma, "higher-order terms:\n");
        for (k = 4; k <= best_fit.order_of_fit; k++)
          fprintf(fpo_ma, "        %.16le +/-  %.16le\n",
                  best_fit.coefficient[k],
                  best_fit.coefficient_error[k]);
      }
      if (verbose > 1)
        fflush(fpo_ma);
    }
  }
  for (i = 0; i < 6; i++)
    if (n_fits[i])
      C[i] /= n_fits[i];
  if (best_fit.coefficient) {
    free(best_fit.coefficient);
  }
  if (best_fit.coefficient_error) {
    free(best_fit.coefficient_error);
  }
  best_fit.coefficient = best_fit.coefficient_error = NULL;

#if DEBUG1
  checkAssignedErrors(fpo_ma, initial, final, error, n_points_total);
#endif

  /* Subtract terms due to Ci, Rij, Tijj, and Qijjj from all final vectors.
   * Also subtract terms due to higher orders in the fits.
   */
  for (k = n_points_total - 1; k >= 0; k--) {
    for (i = 5; i >= 0; i--) {
      final[k][i] -= C[i];
      for (j = 5; j >= 0; j--) {
        if (saved_fit[i][j].order_of_fit > 3) {
          saved_fit[i][j].coefficient[0] = 0;
          if ((term = initial[k][j]))
            final[k][i] -= poly(saved_fit[i][j].coefficient,
                                saved_fit[i][j].order_of_fit + 1, term);
        } else {
          if ((term = initial[k][j]))
            final[k][i] -= term * (R[i][j] + term * (T[i][j][j] + term * Q[i][j][j][j]));
        }
      }
    }
  }
  for (i = 5; i >= 0; i--)
    for (j = 5; j >= 0; j--) {
      free(saved_fit[i][j].coefficient);
      saved_fit[i][j].coefficient = NULL;
      free(saved_fit[i][j].coefficient_error);
      saved_fit[i][j].coefficient_error = NULL;
    }
  free_zarray_2d((void **)saved_fit, 6, 6);

#if DEBUG1
  checkAssignedErrors(fpo_ma, initial, final, error, n_points_total);
#endif

  /* Now compute Tijk's, Qijkk's, and Qijjk's for j>k. */
  for (j = 0; j < 6; j++) {
    if (maximum_value[j] <= 0)
      continue;
    for (k = 0; k < j; k++) {
      if (maximum_value[k] <= resolution[k]) {
        fprintf(fpo_ma,
                "too few pairs for computing T[i][j][%ld] and Q[i][j][%ld][%ld]\n",
                k + 1, k + 1, k + 1);
        continue;
      }
      /* Choose a set of pairs of initial and final vectors such
       * that only j and k coordinates are ever non-zero, and such
       * that j coordinates take their maximum value.
       */
      for (i = 5; i >= 0; i--)
        if (i == k) {
          constraint_code[i] = UNCONSTRAINED;
        } else if (i == j) {
          constraint[i] = maximum_value[i];
          constraint_code[i] = CONSTRAIN_EQ;
        } else {
          constraint[i] = 0.0;
          constraint_code[i] = CONSTRAIN_EQ;
        }

      if ((n_pairs1 = selectInitialVectors(set1_i, set1_f, set1_error,
                                           n_points_max, initial, final, error,
                                           n_points_total, resolution, constraint,
                                           constraint_code)) <= 3) {
        fprintf(fpo_ma,
                "too few pairs for computing T[i][%ld][%ld] and Q[i][%ld][%ld][%ld]\n",
                j + 1, k + 1, j + 1, k + 1, k + 1);
        continue;
      }
      if (verbose > 1)
        fprintf(fpo_ma, "%ld pairs selected for T[i][%ld][%ld] and Q[i][%ld][%ld][%ld]\n",
                n_pairs1, j + 1, k + 1, j + 1, k + 1, k + 1);
#if DEBUG2
      showSelectedPairs(fpo_ma, set1_i, set1_f, set1_error, n_pairs1);
      checkAssignedErrors(fpo_ma, set1_i, set1_f, set1_error, n_pairs1);
#endif

      constraint[j] = -maximum_value[j];
      constraint_code[j] = CONSTRAIN_EQ;

      if ((n_pairs2 = selectInitialVectors(set2_i, set2_f, set2_error,
                                           n_points_max, initial, final, error,
                                           n_points_total, resolution, constraint,
                                           constraint_code)) <= 3) {
        fprintf(fpo_ma,
                "too few pairs for computing T[i][%ld][%ld] and Q[i][%ld][%ld][%ld]\n",
                j + 1, k + 1, j + 1, k + 1, k + 1);
        continue;
      }
      if (verbose > 1)
        fprintf(fpo_ma, "%ld additional pairs selected for T[i][%ld][%ld] and Q[i][%ld][%ld][%ld]\n",
                n_pairs2, j + 1, k + 1, j + 1, k + 1, k + 1);
#if DEBUG2
      fprintf(fpo_ma, "vectors selected:\n");
      showSelectedPairs(fpo_ma, set1_i, set1_f, set1_error, n_pairs1);
      checkAssignedErrors(fpo_ma, initial, final, error, n_points_total);
      checkAssignedErrors(fpo_ma, set2_i, set2_f, set2_error, n_pairs2);
#endif

      for (i = 5; i >= 0; i--) {
        if (!findBestFit(&best_fit1, set1_i, set1_f, set1_error,
                         n_pairs1, i, k, 10 * (10 * i + k) + 2, max_order_of_fits)) {
          fprintf(fpo_ma, "unable to find fit1 for R[%ld][%ld], T[%ld][%ld][%ld], and Q[%ld][%ld][%ld][%ld]\n",
                  i + 1, j + 1, i + 1, j + 1, j + 1,
                  i + 1, j + 1, j + 1, j + 1);
          fflush(fpo_ma);
          continue;
        }

#if DEBUG2
        checkAssignedErrors(fpo_ma, initial, final, error, n_points_total);
        checkAssignedErrors(fpo_ma, set2_i, set2_f, set2_error, n_pairs2);
        fprintf(fpo_ma, "order of best fit1 = %ld\n",
                best_fit1.order_of_fit);
        fprintf(fpo_ma, "chi-squared = %le \n", best_fit1.chi_squared);
        fprintf(fpo_ma, "coefficients of fit1:\n");
        for (l = 0; l < best_fit1.order_of_fit + 1; l++)
          fprintf(fpo_ma, "a[%ld] = %le +/- %le\n",
                  l, best_fit1.coefficient[l], best_fit1.coefficient_error[l]);
        fflush(fpo_ma);
        checkAssignedErrors(fpo_ma, initial, final, error, n_points_total);
        checkAssignedErrors(fpo_ma, set2_i, set2_f, set2_error, n_pairs2);
#endif
        if (!findBestFit(&best_fit2, set2_i, set2_f, set2_error,
                         n_pairs2, i, k, 10 * (10 * i + k) + 3, max_order_of_fits)) {
          fprintf(fpo_ma, "unable to find fit2 for R[%ld][%ld], T[%ld][%ld][%ld], and Q[%ld][%ld][%ld][%ld]\n",
                  i + 1, j + 1, i + 1, j + 1, j + 1,
                  i + 1, j + 1, j + 1, j + 1);
          fflush(fpo_ma);
          continue;
        }

        if (verbose) {
          fprintf(fpo_ma, "order of fits for for computing T[%ld][%ld][%ld] and Q[%ld][%ld][%ld][%ld] = %d, %d\n",
                  i + 1, j + 1, k + 1, i + 1, j + 1, k + 1, k + 1, best_fit1.order_of_fit, best_fit2.order_of_fit);
        }
#if DEBUG2
        fprintf(fpo_ma, "order of best fit2 = %d\n",
                best_fit2.order_of_fit);
        fprintf(fpo_ma, "chi-squared = %le \n", best_fit2.chi_squared);
        fprintf(fpo_ma, "coefficients of fit2:\n");
        for (l = 0; l < best_fit2.order_of_fit + 1; l++)
          fprintf(fpo_ma, "a[%ld] = %le +/- %le\n",
                  l, best_fit2.coefficient[l], best_fit2.coefficient_error[l]);
        fflush(fpo_ma);
        checkAssignedErrors(fpo_ma, initial, final, error, n_points_total);
#endif
        T[i][j][k] =
          (best_fit1.coefficient[1] - best_fit2.coefficient[1]) / maximum_value[j] / 2.0;
        Q[i][j][j][k] =
          (best_fit1.coefficient[1] + best_fit2.coefficient[1]) / sqr(maximum_value[j]) / 2.0;
        if (verbose > 1) {
          fprintf(fpo_ma, "    T[%ld][%ld][%ld] = %.16le +/- %.16le\n",
                  i + 1, j + 1, k + 1, T[i][j][k],
                  best_fit1.coefficient_error[1] / maximum_value[j] / sqrt(2.0));
          fprintf(fpo_ma, "    Q[%ld][%ld][%ld][%ld] = %.16le +/- %.16le\n",
                  i + 1, j + 1, j + 1, k + 1, Q[i][j][j][k],
                  best_fit1.coefficient_error[1] / sqr(maximum_value[j]) / sqrt(2.0));
        }
        if (best_fit1.order_of_fit >= 2 && best_fit2.order_of_fit >= 2) {
          Q[i][j][k][k] =
            (best_fit1.coefficient[2] - best_fit2.coefficient[2]) / maximum_value[j] / 2.0;
          if (verbose)
            fprintf(fpo_ma, "    Q[%ld][%ld][%ld][%ld] = %.16le +/- %.16le\n",
                    i + 1, j + 1, k + 1, k + 1, Q[i][j][k][k],
                    best_fit1.coefficient_error[2] / maximum_value[j] / sqrt(2.0));
        }
        fflush(fpo_ma);
      }
      /* compute correction to Tijk from fourth order term Vijkk */
#if DEBUG2
      if (verbose)
        fprintf(fpo_ma, "computing corrections to Tijk's\n");
#endif
      /* Choose a set of pairs of initial and final vectors such
       * that only j and k coordinates are ever non-zero, and such
       * that j coordinates take their nearest-to-maximum value.
       */
      for (i = 5; i >= 0; i--)
        if (i == k) {
          constraint_code[i] = UNCONSTRAINED;
        } else if (i == j) {
          constraint[i] = maximum_value[i] - step_size[i];
          constraint_code[i] = CONSTRAIN_EQ;
        } else {
          constraint[i] = 0.0;
          constraint_code[i] = CONSTRAIN_EQ;
        }

      if ((n_pairs1 = selectInitialVectors(set1_i, set1_f, set1_error,
                                           n_points_max, initial, final, error,
                                           n_points_total, resolution, constraint,
                                           constraint_code)) <= 3) {
        fprintf(fpo_ma,
                "too few pairs for computing T[i][%ld][%ld] and Q[i][%ld][%ld][%ld]\n",
                j + 1, k + 1, j + 1, k + 1, k + 1);
        continue;
      }
      if (verbose)
        fprintf(fpo_ma, "%ld pairs selected for T[i][%ld][%ld] and Q[i][%ld][%ld][%ld]\n",
                n_pairs1, j + 1, k + 1, j + 1, k + 1, k + 1);
#if DEBUG2
      showSelectedPairs(fpo_ma, set1_i, set1_f, set1_error, n_pairs1);
#endif

      constraint[j] = -constraint[j];
      constraint_code[j] = CONSTRAIN_EQ;

      if ((n_pairs2 = selectInitialVectors(set2_i, set2_f, set2_error,
                                           n_points_max, initial, final, error,
                                           n_points_total, resolution, constraint,
                                           constraint_code)) <= 3) {
        fprintf(fpo_ma,
                "too few pairs for computing T[i][%ld][%ld] and Q[i][%ld][%ld][%ld]\n",
                j + 1, k + 1, j + 1, k + 1, k + 1);
        fflush(fpo_ma);
        continue;
      }
      if (verbose)
        fprintf(fpo_ma, "%ld additional pairs selected for T[i][%ld][%ld] and Q[i][%ld][%ld][%ld]\n",
                n_pairs2, j + 1, k + 1, j + 1, k + 1, k + 1);
#if DEBUG2
      fprintf(fpo_ma, "vectors selected:\n");
      showSelectedPairs(fpo_ma, set1_i, set1_f, set1_error, n_pairs1);
      checkAssignedErrors(fpo_ma, initial, final, error, n_points_total);
#endif

      for (i = 5; i >= 0; i--) {
        if (!findBestFit(&best_fit1, set1_i, set1_f, set1_error,
                         n_pairs1, i, k, 10 * (10 * i + k) + 2, max_order_of_fits)) {
          fprintf(fpo_ma, "unable to find fit1 for R[%ld][%ld], T[%ld][%ld][%ld], and Q[%ld][%ld][%ld][%ld]\n",
                  i + 1, j + 1, i + 1, j + 1, j + 1,
                  i + 1, j + 1, j + 1, j + 1);
          fflush(fpo_ma);
          continue;
        }

#if DEBUG2
        fprintf(fpo_ma, "order of best fit1 = %ld\n",
                best_fit1.order_of_fit);
        fprintf(fpo_ma, "chi-squared = %le \n", best_fit1.chi_squared);
        fprintf(fpo_ma, "coefficients of fit1:\n");
        for (l = 0; l < best_fit1.order_of_fit + 1; l++)
          fprintf(fpo_ma, "a[%ld] = %le +/- %le\n",
                  l, best_fit1.coefficient[l], best_fit1.coefficient_error[l]);
        fflush(fpo_ma);
        checkAssignedErrors(fpo_ma, initial, final, error, n_points_total);
#endif
        if (!findBestFit(&best_fit2, set2_i, set2_f, set2_error,
                         n_pairs2, i, k, 10 * (10 * i + k) + 3, max_order_of_fits)) {
          fprintf(fpo_ma, "unable to find fit2 for R[%ld][%ld], T[%ld][%ld][%ld], and Q[%ld][%ld][%ld][%ld]\n",
                  i + 1, j + 1, i + 1, j + 1, j + 1,
                  i + 1, j + 1, j + 1, j + 1);
          fflush(fpo_ma);
          continue;
        }
#if DEBUG2
        fprintf(fpo_ma, "order of best fit2 = %ld\n",
                best_fit2.order_of_fit);
        fprintf(fpo_ma, "chi-squared = %le \n", best_fit2.chi_squared);
        fprintf(fpo_ma, "coefficients of fit2:\n");
        for (l = 0; l < best_fit2.order_of_fit + 1; l++)
          fprintf(fpo_ma, "a[%ld] = %le +/- %le\n",
                  l, best_fit2.coefficient[l], best_fit2.coefficient_error[l]);
        fflush(fpo_ma);
        checkAssignedErrors(fpo_ma, initial, final, error, n_points_total);
#endif
        if (verbose > 1)
          fprintf(fpo_ma, "order of fits for for computing T[%ld][%ld][%ld] and Q[%ld][%ld][%ld][%ld] = %d, %d\n",
                  i + 1, j + 1, k + 1, i + 1, j + 1, k + 1, k + 1, best_fit1.order_of_fit, best_fit2.order_of_fit);
        /* compute second value of T[i][j][k] */
        Tijk = (best_fit1.coefficient[1] - best_fit2.coefficient[1]) / (maximum_value[j] - step_size[j]) / 2.0;
        /* compute correction to stored value */
        T[i][j][k] -= sqr(maximum_value[j]) * (T[i][j][k] - Tijk) /
                      ((2 * maximum_value[j] - step_size[j]) * step_size[j]);
        if (verbose > 1) {
          fprintf(fpo_ma, "T[%ld][%ld][%ld] = %.16le    <-- second approx\n", i + 1, j + 1, k + 1, Tijk);
          fprintf(fpo_ma, "T[%ld][%ld][%ld] = %.16le    <-- corrected\n",
                  i + 1, j + 1, k + 1, T[i][j][k]);
          fflush(fpo_ma);
        }
      }
    }
  }
  if (best_fit1.coefficient) {
    free(best_fit1.coefficient);
    best_fit1.coefficient = NULL;
  }
  if (best_fit2.coefficient) {
    free(best_fit2.coefficient);
    best_fit2.coefficient = NULL;
  }
  if (best_fit1.coefficient_error) {
    free(best_fit1.coefficient_error);
    best_fit1.coefficient_error = NULL;
  }
  if (best_fit2.coefficient_error) {
    free(best_fit2.coefficient_error);
    best_fit2.coefficient_error = NULL;
  }

  /* Subtract Tijk, Qijkk and Qijjk terms (for j>k) from final vectors. */
  for (i = 5; i >= 0; i--) {
    for (l = n_points_total - 1; l >= 0; l--) {
      term = final[l][i];
      for (j = 0; j < 6; j++)
        for (k = 0; k < j; k++)
          term -= T[i][j][k] * initial[l][j] * initial[l][k] + Q[i][j][k][k] * initial[l][j] * sqr(initial[l][k]) + Q[i][j][j][k] * sqr(initial[l][j]) * initial[l][k];
      final[l][i] = term;
    }
  }

  /* Now compute Qijkl terms for j!=k, k!=l, l!=j */
  if (2 > n_points_max) {
    free_zarray_2d((void **)set1_i, n_points_max, 6);
    free_zarray_2d((void **)set1_f, n_points_max, 6);
    free_zarray_2d((void **)set1_error, n_points_max, 6);
    free_zarray_2d((void **)set2_i, n_points_max, 6);
    free_zarray_2d((void **)set2_f, n_points_max, 6);
    free_zarray_2d((void **)set2_error, n_points_max, 6);
    n_points_max = 2;
    set1_i = (double **)zarray_2d(sizeof(**set1_i), n_points_max, 6);
    set1_f = (double **)zarray_2d(sizeof(**set1_f), n_points_max, 6);
    set1_error = (double **)zarray_2d(sizeof(**set1_error), n_points_max, 6);
    set2_i = (double **)zarray_2d(sizeof(**set2_i), n_points_max, 6);
    set2_f = (double **)zarray_2d(sizeof(**set2_f), n_points_max, 6);
    set2_error = (double **)zarray_2d(sizeof(**set2_error), n_points_max, 6);
  }

  for (j = 0; j < 6; j++) {
    if (maximum_value[j] == 0)
      continue;
    for (k = 0; k < j; k++) {
      if (maximum_value[k] == 0)
        continue;
      for (l = 0; l < k; l++) {
        if (maximum_value[l] == 0)
          continue;
        /* Choose a set of pairs of initial and final vectors such
	 * that only j, k, and l coordinates are non-zero, and each
	 * non-zero component takes its maximum value
	 */
        for (i = 5; i >= 0; i--)
          if (i != j && i != k && i != l) {
            constraint[i] = 0.0;
            constraint_code[i] = CONSTRAIN_EQ;
          } else {
            constraint[i] = maximum_value[i];
            constraint_code[i] = CONSTRAIN_EQ;
          }

        if ((n_pairs1 = selectInitialVectors(set1_i, set1_f, set1_error,
                                             n_points_max, initial, final, error,
                                             n_points_total, resolution, constraint,
                                             constraint_code)) != 1) {
          fprintf(fpo_ma,
                  "wrong number of pairs for computing Q[i][%ld][%ld][%ld]\n",
                  j + 1, k + 1, l + 1);
          continue;
        }
#if DEBUG3
        fprintf(fpo_ma, "first vectors selected for j=%ld, k=%ld, l=%ld\n",
                j, k, l);
        for (i = 0; i < n_pairs1; i++) {
          fprintf(fpo_ma, "ini: ");
          for (m = 0; m < 6; m++)
            fprintf(fpo_ma, "%.16le ", set1_i[i][m]);
          fprintf(fpo_ma, "\nfin: ");
          for (m = 0; m < 6; m++)
            fprintf(fpo_ma, "%.16le ", set1_f[i][m]);
          fprintf(fpo_ma, "\n");
        }
        fflush(fpo_ma);
#endif
        /* Choose a set of pairs of initial and final vectors such
	 * that only j, k, and l coordinates are non-zero, and each
	 * non-zero component takes its minimum value
	 */
        for (i = 5; i >= 0; i--)
          if (i != j && i != k && i != l) {
            constraint[i] = 0.0;
            constraint_code[i] = CONSTRAIN_EQ;
          } else {
            constraint[i] = -maximum_value[i];
            constraint_code[i] = CONSTRAIN_EQ;
          }

        if ((n_pairs2 = selectInitialVectors(set2_i, set2_f, set2_error,
                                             n_points_max, initial, final, error,
                                             n_points_total, resolution, constraint,
                                             constraint_code)) != 1) {
          fprintf(fpo_ma,
                  "wrong number of pairs for computing Q[i][%ld][%ld][%ld]\n",
                  j + 1, k + 1, l + 1);
          continue;
        }
#if DEBUG3
        fprintf(fpo_ma, "second vectors selected for j=%ld, k=%ld, l=%ld\n",
                j, k, l);
        for (i = 0; i < n_pairs1; i++) {
          fprintf(fpo_ma, "ini: ");
          for (m = 0; m < 6; m++)
            fprintf(fpo_ma, "%.16le ", set2_i[i][m]);
          fprintf(fpo_ma, "\nfin: ");
          for (m = 0; m < 6; m++)
            fprintf(fpo_ma, "%.16le ", set2_f[i][m]);
          fprintf(fpo_ma, "\n");
        }
        fflush(fpo_ma);
#endif
        for (i = 5; i >= 0; i--) {
          Q[i][j][k][l] = (set1_f[0][i] - set2_f[0][i]) / (2 * set1_i[0][j] * set1_i[0][k] * set1_i[0][l]);
          if (verbose > 1)
            fprintf(fpo_ma, "Q[%ld][%ld][%ld][%ld] = %le\n",
                    i + 1, j + 1, k + 1, l + 1, Q[i][j][k][l]);
        }
      }
    }
  }

  /* Subtract Qijkl terms (for j>k, k>l ) from final vectors. */
  for (i = 0; i < 6; i++) {
    sum = 0;
    for (m = n_points_total - 1; m >= 0; m--) {
      term = final[m][i];
      for (j = 0; j < 6; j++) {
        for (k = 0; k < j; k++) {
          for (l = 0; l < k; l++)
            term -= Q[i][j][k][l] * initial[m][j] * initial[m][k] * initial[m][l];
        }
      }
      sum += fabs(term);
    }
    if (verbose > 1)
      fprintf(fpo_ma, "average absolute residual of fit for %ldth coordinate: %le\n",
              i, sum / n_points_total);
  }

  free(resolution);
  free(constraint);
  free(n_fits);
  free(constraint_code);
  free_zarray_2d((void **)set1_i, n_points_max, 6);
  free_zarray_2d((void **)set1_f, n_points_max, 6);
  free_zarray_2d((void **)set1_error, n_points_max, 6);
  free_zarray_2d((void **)set2_i, n_points_max, 6);
  free_zarray_2d((void **)set2_f, n_points_max, 6);
  free_zarray_2d((void **)set2_error, n_points_max, 6);

  return M;
}

int selectInitialVectors(
  double **set_i, double **set_f, double **set_e,
  int n_max_select, double **initial, double **final, double **error,
  int n_vectors, double *resolution, double *constraint, int *constraint_code) {
  register int i_component, acceptable;
  register double *vector;
  int i_pair, n_pairs;

  n_pairs = 0;
  for (i_pair = n_vectors - 1; n_pairs < n_max_select && i_pair >= 0; i_pair--) {
    acceptable = 1;
    vector = initial[i_pair];
    for (i_component = 5; acceptable && i_component >= 0; i_component--) {
      switch (constraint_code[i_component]) {
      case UNCONSTRAINED:
        break;
      case CONSTRAIN_EQ:
        if (fabs(vector[i_component] - constraint[i_component]) > resolution[i_component])
          acceptable = 0;
        break;
      case CONSTRAIN_NEQ:
        if (fabs(vector[i_component] - constraint[i_component]) <= resolution[i_component])
          acceptable = 0;
        break;
      default:
        printf("unknown constraint code %d\n", constraint_code[i_component]);
        exit(1);
        break;
      }
    }
    if (acceptable) {
      for (i_component = 5; i_component >= 0; i_component--) {
        set_i[n_pairs][i_component] = vector[i_component];
        set_f[n_pairs][i_component] = final[i_pair][i_component];
        set_e[n_pairs][i_component] = error[i_pair][i_component];
      }
      n_pairs++;
    }
  }
  return (n_pairs);
}

int findBestFit(FIT *best_fit, double **set_i, double **set_f, double **set_error,
                int n_pairs, int i, int j, /* as in Tijk and Qijkl */
                int type, int max_order) {
  FIT fit;
  register int k, order_of_fit, lowest_order_of_fit, highest_order_of_fit;
  double *xdata, *ydata, *sigmay;
  long nErrorsZero;

  if (n_pairs < 3)
    return 0;
  lowest_order_of_fit = 2;
  highest_order_of_fit = (max_order > n_pairs - 2 ? n_pairs - 2 : max_order);
  xdata = tmalloc(sizeof(*xdata) * n_pairs);
  ydata = tmalloc(sizeof(*ydata) * n_pairs);
  sigmay = tmalloc(sizeof(*sigmay) * n_pairs);

  fit.coefficient = tmalloc((2 * highest_order_of_fit + 1) * sizeof(*(fit.coefficient)));
  fit.coefficient_error = tmalloc((2 * highest_order_of_fit + 1) * sizeof(*(fit.coefficient_error)));

  if (best_fit->coefficient) {
    free(best_fit->coefficient);
    best_fit->coefficient = NULL;
  }
  if (best_fit->coefficient_error) {
    free(best_fit->coefficient_error);
    best_fit->coefficient_error = NULL;
  }
  best_fit->order_of_fit = 0;

  nErrorsZero = 0;
  for (k = n_pairs - 1; k >= 0; k--) {
    xdata[k] = set_i[k][j];
    ydata[k] = set_f[k][i];
    if ((sigmay[k] = set_error[k][i]) == 0) {
      /*
      printf("sigma[%d] = 0 in findBestFit()\n", k);
      showSelectedPairs(stdout, set_i, set_f, set_error, n_pairs);
      */
      nErrorsZero++;
    }
  }
  if (nErrorsZero) {
    if (nErrorsZero != n_pairs) {
      printf("Fatal error: some but not all integration errors are zero in findBestFit()\n");
      exit(1);
    }
    /* Just do an equal-weights fit */
    for (k = n_pairs - 1; k >= 0; k--)
      sigmay[k] = 1;
  }

  best_fit->chi_squared = HUGE_VAL;
  best_fit->order_of_fit = -1;
  order_of_fit = lowest_order_of_fit;

  do {
    if (!lsfn(xdata, ydata, sigmay, n_pairs, order_of_fit,
              fit.coefficient, fit.coefficient_error, &(fit.chi_squared), NULL)) {
      printf("fitting error: order %d, type %d\n", order_of_fit, type);
      continue;
    }
    if (fit.chi_squared < best_fit->chi_squared) {
      best_fit->chi_squared = fit.chi_squared;
      if (best_fit->coefficient != NULL) {
        free(best_fit->coefficient);
        best_fit->coefficient = NULL;
        free(best_fit->coefficient_error);
        best_fit->coefficient_error = NULL;
      }
      best_fit->coefficient = fit.coefficient;
      best_fit->coefficient_error = fit.coefficient_error;
      fit.coefficient = NULL;
      fit.coefficient_error = NULL;
      best_fit->order_of_fit = order_of_fit;
      fit.coefficient = tmalloc((2 * highest_order_of_fit + 1) * sizeof(*(fit.coefficient)));
      fit.coefficient_error = tmalloc((2 * highest_order_of_fit + 1) * sizeof(*(fit.coefficient_error)));
    }
  } while (++order_of_fit <= highest_order_of_fit);

  free(xdata);
  free(ydata);
  free(sigmay);
  free(fit.coefficient);
  fit.coefficient = NULL;
  free(fit.coefficient_error);
  fit.coefficient_error = NULL;

  if (best_fit->order_of_fit == -1)
    return 0;

  for (k = 0; k <= best_fit->order_of_fit; k++) {
    if (fabs(best_fit->coefficient[k]) < best_fit->coefficient_error[k])
      best_fit->coefficient[k] = 0;
  }

  return (best_fit->order_of_fit);
}

void copyFit(FIT *target, FIT *source) {
  int i;

  if (target->coefficient)
    free(target->coefficient);
  target->coefficient = array_1d(sizeof(*target->coefficient), 0, source->order_of_fit + 1);
  if (target->coefficient_error)
    free(target->coefficient_error);
  target->coefficient_error = array_1d(sizeof(*target->coefficient_error), 0, source->order_of_fit + 1);
  target->order_of_fit = source->order_of_fit;
  target->chi_squared = source->chi_squared;

  for (i = 0; i <= source->order_of_fit; i++) {
    target->coefficient[i] = source->coefficient[i];
    target->coefficient_error[i] = source->coefficient_error[i];
  }
}

void showSelectedPairs(FILE *fp, double **ini, double **fin, double **err, int n) {
  int i;
  for (i = 0; i < n; i++) {
    fprintf(fp, "%3d ini: %le %le %le %le %le %le\n", i,
            ini[i][0], ini[i][1], ini[i][2], ini[i][3], ini[i][4], ini[i][5]);
    fprintf(fp, "    fin: %le %le %le %le %le %le\n",
            fin[i][0], fin[i][1], fin[i][2], fin[i][3], fin[i][4], fin[i][5]);
    fprintf(fp, "    err: %le %le %le %le %le %le\n",
            err[i][0], err[i][1], err[i][2], err[i][3], err[i][4], err[i][5]);
  }
}

void checkAssignedErrors(FILE *fp, double **ini, double **fin, double **err, int n) {
  int i, j;
  for (i = 0; i < n; i++) {
    for (j = 0; j < 6; j++) {
      if (err[i][j] != err[0][0]) {
        fprintf(fp, "bad error value for vector %d, component %d\n", i, j);
        showSelectedPairs(fp, ini, fin, err, n);
        return;
      }
    }
  }
  fprintf(fp, "error vectors are okay\n");
}
