/*************************************************************************\
* Copyright (c) 2019 The University of Chicago, as Operator of Argonne
* National Laboratory.
* Copyright (c) 2019 The Regents of the University of California, as
* Operator of Los Alamos National Laboratory.
* This file is distributed subject to a Software License Agreement found
* in the file LICENSE that is included with this distribution. 
\*************************************************************************/

/*
 * Michael Borland, 2019
 */
#include "mdb.h"
#include "track.h"

double trimApertureSearchResult1(long lines,
                                 double *xLimit, double *yLimit,
                                 double *dxFactor, double *dyFactor,
                                 long fullPlane) {
  long line;
  double area;
  if (!fullPlane) {
    /* First clip off any portions that stick out like islands.  This is done in three steps. */

    /* 1.  Insist that the x values must be monotonically increasing 
     */
    for (line = 0; line < lines / 2; line++)
      if (xLimit[line + 1] < xLimit[line])
        xLimit[line + 1] = xLimit[line];
    for (line = lines - 1; line > lines / 2; line--)
      if (xLimit[line - 1] > xLimit[line])
        xLimit[line - 1] = xLimit[line];

    /* 2. for x<0, y values must increase monotonically as x increases (larger index) */
    for (line = lines / 2; line > 0; line--) {
      if (yLimit[line - 1] > yLimit[line])
        yLimit[line - 1] = yLimit[line];
    }

    /* 3. for x>0, y values must fall monotonically as x increases (larger index) */
    for (line = lines / 2; line < (lines - 1); line++) {
      if (yLimit[line + 1] > yLimit[line])
        yLimit[line + 1] = yLimit[line];
    }

  } else { /* full plane */
    /* First clip off any portions that stick out like islands.  */

    for (line = 0; line < lines - 1; line++) {
      if (dxFactor[line] > 0 && dxFactor[line + 1] > 0) {
        if (dyFactor[line] > 0 && dyFactor[line + 1] > 0) {
          /* Quadrant 1 */
          if (xLimit[line + 1] < xLimit[line])
            xLimit[line] = xLimit[line + 1];
        } else if (dyFactor[line] < 0 && dyFactor[line + 1] < 0) {
          /* Quadrant 4 */
          if (xLimit[line + 1] > xLimit[line])
            xLimit[line + 1] = xLimit[line];
        }
      } else if (dxFactor[line] < 0 && dxFactor[line + 1] < 0) {
        if (dyFactor[line] > 0 && dyFactor[line + 1] > 0) {
          /* Quadrant 2 */
          if (xLimit[line + 1] < xLimit[line])
            xLimit[line + 1] = xLimit[line];
        } else if (dyFactor[line] < 0 && dyFactor[line + 1] < 0) {
          /* Quadrant 3 */
          if (xLimit[line + 1] > xLimit[line])
            xLimit[line] = xLimit[line + 1];
        }
      }
    }

    for (line = 0; line < lines - 1; line++) {
      if (dxFactor[line] > 0 && dxFactor[line + 1] > 0) {
        if (dyFactor[line] > 0 && dyFactor[line + 1] > 0) {
          /* Quadrant 1 */
          if (yLimit[line + 1] > yLimit[line])
            yLimit[line + 1] = yLimit[line];
        } else if (dyFactor[line] < 0 && dyFactor[line + 1] < 0) {
          /* Quadrant 4 */
          if (yLimit[line] < yLimit[line + 1])
            yLimit[line] = yLimit[line + 1];
        }
      } else if (dxFactor[line] < 0 && dxFactor[line + 1] < 0) {
        if (dyFactor[line] > 0 && dyFactor[line + 1] > 0) {
          /* Quadrant 2 */
          if (yLimit[line + 1] < yLimit[line])
            yLimit[line] = yLimit[line + 1];
        } else if (dyFactor[line] < 0 && dyFactor[line + 1] < 0) {
          /* Quadrant 3 */
          if (yLimit[line + 1] < yLimit[line])
            yLimit[line + 1] = yLimit[line];
        }
      }
    }
  }

  /* perform trapazoid rule integration */
  area = 0;
  for (line = 0; line < lines - 1; line++)
    area += (xLimit[line + 1] - xLimit[line]) * (yLimit[line + 1] + yLimit[line]) / 2;

  return area;
}

double trimApertureSearchResult(long lines,
                                double *xLimit, double *yLimit,
                                double *dxFactor, double *dyFactor,
                                long fullPlane) {
  double area1 = -1, area2 = -1;
  long count;
  count = lines;
  while (count-- &&
         (area2 = trimApertureSearchResult1(lines, xLimit, yLimit, dxFactor, dyFactor, fullPlane)) != area1)
    area1 = area2;
  return area1;
}
