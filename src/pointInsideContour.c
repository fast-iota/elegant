/* This code found on multiple internet sources, translated to C */

#include <stdio.h>
#include <math.h>

double pointIsLeftOfLine(long i1, long i2,     /* indices of vertex points to test */
                         double *x, double *y, /* (x, y) of vertices */
                         double x0, double y0  /* (x, y) of point to test */
);

int pointIsInsideContour(
  double x0, double y0, /* point to check */
  double *x, double *y, /* points on contour */
  long n,               /* points in contour */
  double *center,       /* center of rotation of contour */
  double theta          /* clockwise angle to rotation of contour */
) {
  long winding_number = 0;
  long i1, i2;

  if (theta != 0 && center) {
    double x1, y1, ct, st;
    /* rotate (x0, y0) by theta about center */
    x0 -= center[0];
    y0 -= center[1];
    ct = cos(theta);
    st = sin(theta);
    x1 = x0 * ct - y0 * st + center[0];
    y1 = x0 * st + y0 * ct + center[1];
    x0 = x1;
    y0 = y1;
  }

  for (i1 = 0; i1 < n; i1++) {
    if (i1 == (n - 1))
      /* wrap */
      i2 = 0;
    else
      i2 = i1 + 1;

    if (y[i1] <= y0) {
      /* start y <= y0 */
      if (y[i2] > y0) {
        /* upward crossing */
#ifdef DEBUG
        printf("upward crossing\n");
        fflush(stdout);
#endif
        if (pointIsLeftOfLine(i1, i2, x, y, x0, y0) > 0) {
          /* Point left of edge */
          ++winding_number;
        }
      }
    } else {
      /* start y > y0 */
      if (y[i2] <= y0) {
        /* downward crossing */
#ifdef DEBUG
        printf("downward crossing\n");
        fflush(stdout);
#endif
        if (pointIsLeftOfLine(i1, i2, x, y, x0, y0) < 0) {
          /* Point right of edge */
          --winding_number;
        }
      }
    }
#ifdef DEBUG
    printf("i1 = %ld, winding_number = %ld\n", i1, winding_number);
    fflush(stdout);
#endif
  }

#ifdef DEBUG
  printf("winding_number = %ld\n", winding_number);
  fflush(stdout);
#endif

  return (winding_number != 0);
}

double pointIsLeftOfLine(long i1, long i2,     /* indices of vertex points to test */
                         double *x, double *y, /* (x, y) of vertices */
                         double x0, double y0  /* (x, y) of point to test */
) {
  double result;
  result = (x[i2] - x[i1]) * (y0 - y[i1]) - (x0 - x[i1]) * (y[i2] - y[i1]);
#ifdef DEBUG
  printf("pointIsLeftOfLine(%ld, %ld, %le, %le) = %le\n",
         i1, i2, x0, y0, result);
  fflush(stdout);
#endif

  return result;
}
