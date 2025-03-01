/* Copyright 1994, 2015 by Michael Borland and Argonne National Laboratory,
 * all rights reserved.
 */
/* routine: makeInitialParticleEnsemble()
 * purpose: make ensemble of initial vectors for computation of
 *          transport matrices
 *
 * Michael Borland, 1990.
 */

#include "mdb.h"
#include "track.h"

int makeInitialParticleEnsemble(
  double ***initial, /* will be allocated and initialized */
  double *reference, /* may be provided */
  double ***final,   /* will be allocated and initialized to same values as initial */
  double ***error,   /* will be allocated and initialized */
  int n_points1,     /* number of points in each dimension */
  double *step       /* step size for each dimension */
) {
  int i, j, k, l, i_vector;
  int J, K, L;
  int n_points_total;
  //int n_duplicates;

  /* count up the number of vectors to be integrated */
  /* the unit vectors are e_J, e_K, e_L */
  n_points_total = 1;
  for (J = 5; J >= 0; J--) {
    for (K = J - 1; K >= 0; K--) {
      for (L = K - 1; L >= 0; L--) {
        n_points_total += n_points1 * n_points1 * n_points1 - 1;
      }
    }
  }

  *initial = (double **)czarray_2d(sizeof(***initial), n_points_total, totalPropertiesPerParticle);
  *final = (double **)czarray_2d(sizeof(***final), n_points_total, totalPropertiesPerParticle);
  *error = (double **)czarray_2d(sizeof(***error), n_points_total, totalPropertiesPerParticle);

  /* generate the origin vector */
  for (i = 0; i < 6; i++)
    (*initial)[0][i] = reference ? reference[i] : 0;
  i_vector = 1;

  /* generate other initial vectors */
  for (J = 5; J >= 0; J--) {
    for (K = J - 1; K >= 0; K--) {
      for (L = K - 1; L >= 0; L--) {
        /* J, K, L is the triplet of indices of the unit vectors */
        for (j = (n_points1 - 1) / 2; j >= -(n_points1 - 1) / 2; j--) {
          for (k = (n_points1 - 1) / 2; k >= -(n_points1 - 1) / 2; k--) {
            for (l = (n_points1 - 1) / 2; l >= -(n_points1 - 1) / 2; l--) {
              /* generate initial vector :
	       * j a[J] e[J] + k a[K] e[K] + l a[L] e[L] 
	       */
              if (!j && !k && !l)
                continue; /* don't repeat origin vector */
              if (i_vector >= n_points_total)
                bombElegant("counting problem in make_ensemble()", NULL);
              /* initialize the vector */
              for (i = 5; i >= 0; i--)
                (*initial)[i_vector][i] = reference ? reference[i] : 0;
              /* increment the selected initial coordinates */
              (*initial)[i_vector][J] += j * step[J];
              (*initial)[i_vector][K] += k * step[K];
              (*initial)[i_vector][L] += l * step[L];
              for (i = 5; i >= 0; i--)
                (*final)[i_vector][i] = (*error)[i_vector][i] = 0;
              i_vector++;
            }
          }
        }
      }
    }
  }

  /* check for duplicates */
  //n_duplicates = 0;
  for (i = 0; i < n_points_total; i++) {
    for (j = i + 1; j < n_points_total; j++) {
      if ((*initial)[i][0] == (*initial)[j][0] &&
          (*initial)[i][1] == (*initial)[j][1] &&
          (*initial)[i][2] == (*initial)[j][2] &&
          (*initial)[i][3] == (*initial)[j][3] &&
          (*initial)[i][4] == (*initial)[j][4] &&
          (*initial)[i][5] == (*initial)[j][5]) {
        if (j != (n_points_total - 1))
          memcpy((*initial)[j], (*initial)[n_points_total - 1], totalPropertiesPerParticle * sizeof(double));
        j -= 1;
        n_points_total -= 1;
        //n_duplicates += 1;
      }
    }
  }

  for (i = 0; i < n_points_total; i++) {
    (*initial)[i][6] = i + 1;
    /* since do_tracking() doesn't use initial and final arrays, need to copy the initial coordinates */
    memcpy((*final)[i], (*initial)[i], totalPropertiesPerParticle * sizeof(double));
  }

  /*
  printf("%ld duplicates were removed from initial ensemble, leaving %ld vectors\n", (long)n_duplicates, (long)n_points_total);
  */

  return n_points_total;
}
