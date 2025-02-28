/*************************************************************************\
* Copyright (c) 2002 The University of Chicago, as Operator of Argonne
* National Laboratory.
* Copyright (c) 2002 The Regents of the University of California, as
* Operator of Los Alamos National Laboratory.
* This file is distributed subject to a Software License Agreement found
* in the file LICENSE that is included with this distribution. 
\*************************************************************************/

/* file: chaosMap.c
 * purpose: Do chaos map tracking and analysis.
 *          See file chaosMap.nl for input parameters.
 *
 * Michael Borland, 2004
 */
#include "mdb.h"
#include "track.h"
#include "chaosMap.h"

#define IC_X 0
#define IC_Y 1
#define IC_DELTA 2
#define IC_SURVIVED 3
#define IC_DJX 4
#define IC_DJY 5
#define IC_LOGDJX 6
#define IC_LOGDJY 7
#define N_COLUMNS1 8
static SDDS_DEFINITION column_definition1[N_COLUMNS1] = {
  {"x", "&column name=x, symbol=x, units=m, type=double &end"},
  {"y", "&column name=y, symbol=y, units=m, type=double &end"},
  {"delta", "&column name=delta, type=double &end"},
  {"Survived", "&column name=Survived, type=short &end"},
  {"dJx", "&column name=dJx, units=m, type=double &end"},
  {"dJy", "&column name=dJy, units=m, type=double &end"},
  {"Log10dJx", "&column name=Log10dJx, type=double &end"},
  {"Log10dJy", "&column name=Log10dJy, type=double &end"},
};

#define IC_DF 4
#define IC_LOGDF 5
#define N_COLUMNS2 6
static SDDS_DEFINITION column_definition2[N_COLUMNS2] = {
  {"x", "&column name=x, symbol=x, units=m, type=double &end"},
  {"y", "&column name=y, symbol=y, units=m, type=double &end"},
  {"delta", "&column name=delta, type=double &end"},
  {"Survived", "&column name=Survived, type=short &end"},
  {"dF", "&column name=dF, type=double &end"},
  {"Log10dF", "&column name=Log10dF, type=double &end"},
};

#define IP_STEP 0
#define N_PARAMETERS 2
static SDDS_DEFINITION parameter_definition[N_PARAMETERS] = {
  {"Step", "&parameter name=Step, type=long, description=\"Simulation step\" &end"},
  {"SVNVersion", "&parameter name=SVNVersion, type=string, description=\"SVN version number\", fixed_value=" SVN_VERSION " &end"},
};

static SDDS_DATASET SDDS_cmap;

void setupChaosMap(
  NAMELIST_TEXT *nltext,
  RUN *run,
  VARY *control) {
  /* process namelist input */
  set_namelist_processing_flags(STICKY_NAMELIST_DEFAULTS);
  set_print_namelist_flags(0);
  if (processNamelist(&chaos_map, nltext) == NAMELIST_ERROR)
    bombElegant(NULL, NULL);
  if (echoNamelists)
    print_namelist(stdout, &chaos_map);

  /* check for data errors */
  if (!output)
    bombElegant("no output filename specified", NULL);
  if (xmin > xmax)
    bombElegant("xmin > xmax", NULL);
  if (ymin > ymax)
    bombElegant("ymin > ymax", NULL);
  if (delta_min > delta_max)
    bombElegant("delta_min > delta_max", NULL);
  if (nx < 1)
    nx = 1;
  if (ny < 1)
    ny = 1;
  if (ndelta < 1)
    ndelta = 1;
  if (change_x == 0 || change_y == 0)
    bombElegant("must have change_x and change_y nonzero", NULL);

  output = compose_filename(output, run->rootname);
#if SDDS_MPI_IO
  SDDS_cmap.parallel_io = 1;
  SDDS_MPI_Setup(&SDDS_cmap, 1, n_processors, myid, MPI_COMM_WORLD, 1);
#endif
  if (forward_backward > 0)
    SDDS_ElegantOutputSetup(&SDDS_cmap, output, SDDS_BINARY, 1, "chaos map analysis",
                            run->runfile, run->lattice, parameter_definition, N_PARAMETERS,
                            column_definition2, N_COLUMNS2,
                            "setup_chaosMap", SDDS_EOS_NEWFILE);
  else
    SDDS_ElegantOutputSetup(&SDDS_cmap, output, SDDS_BINARY, 1, "chaos map analysis",
                            run->runfile, run->lattice, parameter_definition, N_PARAMETERS,
                            column_definition1, N_COLUMNS1,
                            "setup_chaosMap", SDDS_EOS_NEWFILE);

  if (control->n_elements_to_vary)
    if (!SDDS_DefineSimpleParameters(&SDDS_cmap, control->n_elements_to_vary,
                                     control->varied_quan_name, control->varied_quan_unit, SDDS_DOUBLE)) {
      SDDS_SetError("Unable to define additional SDDS parameters (setup_chaos_map)");
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors | SDDS_EXIT_PrintErrors);
    }
#if !SDDS_MPI_IO
  if (!SDDS_WriteLayout(&SDDS_cmap)) {
    SDDS_SetError("Unable to write SDDS layout for chaos map");
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors | SDDS_EXIT_PrintErrors);
  }
#endif
}

long doChaosMap(
  RUN *run,
  VARY *control,
  double *referenceCoord,
  ERRORVAL *errcon,
  LINE_LIST *beamline) {
  double startingCoord[6];
  short survived[3];
  double **trackingBuffer = NULL;
  double dx, dy, ddelta, x, y, delta;
  long ix, iy, idelta, ip;
  LINE_LIST *btBeamline; /* back-tracking beamline */
  static double **one_part;
  double p, dJx, dJy, dF;
  long n_part;
#if USE_MPI
  double oldPercentage = 0;
#endif
#if SDDS_MPI_IO
  long points;
  /* Open file here for parallel IO */
  if (!SDDS_MPI_File_Open(SDDS_cmap.MPI_dataset, SDDS_cmap.layout.filename, SDDS_MPI_WRITE_ONLY))
    SDDS_MPI_BOMB("SDDS_MPI_File_Open failed.", &SDDS_cmap.MPI_dataset->MPI_file);
  if (!SDDS_MPI_WriteLayout(&SDDS_cmap))
    SDDS_MPI_BOMB("SDDS_MPI_WriteLayout failed.", &SDDS_cmap.MPI_dataset->MPI_file);

  points = ndelta * nx * ny / n_processors;
  if (myid < (ndelta * nx * ny) % n_processors)
    points++;
  if (!SDDS_StartPage(&SDDS_cmap, points) ||
      !SDDS_SetParameters(&SDDS_cmap, SDDS_SET_BY_INDEX | SDDS_PASS_BY_VALUE, 0, control->i_step, -1)) {
    SDDS_SetError("Unable to start SDDS page (doChaosMap)");
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors | SDDS_EXIT_PrintErrors);
  }
#else
  if (!SDDS_StartPage(&SDDS_cmap, ndelta * nx * ny) ||
      !SDDS_SetParameters(&SDDS_cmap, SDDS_SET_BY_INDEX | SDDS_PASS_BY_VALUE, 0, control->i_step, -1)) {
    SDDS_SetError("Unable to start SDDS page (doChaosMap)");
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors | SDDS_EXIT_PrintErrors);
  }
#endif

  if (!(beamline->flags & BEAMLINE_TWISS_CURRENT))
    bombElegant("Must compute twiss parameters for chaos map\n", NULL);

  btBeamline = NULL;
  if (forward_backward > 0) {
    btBeamline = tmalloc(sizeof(*btBeamline));
    memcpy(btBeamline, beamline, sizeof(*btBeamline));
    copy_line(btBeamline->elem, beamline->elem, beamline->n_elems, 1, NULL, NULL);
    modify_for_backtracking(btBeamline->elem);
  }

#if USE_MPI
  if (verbosity && myid == 1)
    dup2(fdStdout, fileno(stdout)); /* slave will provide warnings etc */
#endif

  if (control->n_elements_to_vary) {
    for (ip = 0; ip < control->n_elements_to_vary; ip++)
      if (!SDDS_SetParameters(&SDDS_cmap, SDDS_SET_BY_INDEX | SDDS_PASS_BY_VALUE, ip + 1,
                              control->varied_quan_value[ip], -1)) {
        SDDS_SetError("Unable to start SDDS page (doChaosMap)");
        SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors | SDDS_EXIT_PrintErrors);
      }
  }

  /* Perform fiducialization by tracking one turn */
  if (!one_part)
    one_part = (double **)czarray_2d(sizeof(**one_part), 1, totalPropertiesPerParticle);
  n_part = 1;
  if (referenceCoord) {
    long i;
    for (i = 0; i < 6; i++)
      one_part[0][i] = referenceCoord[i];
  }
  p = run->p_central;
  if (!do_tracking(NULL, one_part, n_part, NULL, beamline, &p, (double **)NULL, (BEAM_SUMS **)NULL, (long *)NULL,
                   NULL, run, 0, TEST_PARTICLES, 1, 0,
                   NULL, NULL, NULL, NULL, NULL)) {
    printf("Error: lost particle when fiducializing\n");
    exitElegant(1);
  }
  printf("Tracked fiducial particle\n");
  fflush(stdout);

  if (nx > 1)
    dx = (xmax - xmin) / (nx - 1);
  else
    dx = 0;
  if (ny > 1)
    dy = (ymax - ymin) / (ny - 1);
  else
    dy = 0;
  if (ndelta > 1)
    ddelta = (delta_max - delta_min) / (ndelta - 1);
  else
    ddelta = 0;
  ip = 0;
  /* turns = control->n_passes; */

  trackingBuffer = (double **)czarray_2d(sizeof(**trackingBuffer), 3, totalPropertiesPerParticle);

  for (idelta = 0; idelta < ndelta; idelta++) {
    delta = delta_min + idelta * ddelta;
    for (ix = 0; ix < nx; ix++) {
      x = xmin + ix * dx;
      for (iy = 0; iy < ny; iy++) {
        y = ymin + iy * dy;
        memcpy(startingCoord, referenceCoord, sizeof(*startingCoord) * 6);
        startingCoord[0] += x;
        startingCoord[2] += y;
        startingCoord[5] += delta;
#if USE_MPI
        if (myid == (idelta * nx * ny + ix * ny + iy) % n_processors) /* Partition the job according to particle ID */
#endif
        {
          if (forward_backward > 0) {
            int iteration;
            memcpy(trackingBuffer[0], startingCoord, sizeof(*startingCoord) * 6);
            trackingBuffer[0][6] = 1;
            p = run->p_central;
            memset(survived, 0, 3 * sizeof(*survived));
            for (iteration = 0; iteration < forward_backward; iteration++) {
              if (!(survived[0] = do_tracking(NULL, trackingBuffer, 1, NULL, beamline,
                                              &p, (double **)NULL, (BEAM_SUMS **)NULL, (long *)NULL,
                                              NULL, run, 0, TEST_PARTICLES, control->n_passes, 0,
                                              NULL, NULL, NULL, NULL, NULL)))
                break;
              if (!(survived[1] = do_tracking(NULL, trackingBuffer, 1, NULL, btBeamline,
                                              &p, (double **)NULL, (BEAM_SUMS **)NULL, (long *)NULL,
                                              NULL, run, 0, TEST_PARTICLES, control->n_passes, 0,
                                              NULL, NULL, NULL, NULL, NULL)))
                break;
            }
            /* Compute deltas */
            dF = 0;
            if (survived[0] && survived[1]) {
              double beta, alpha, du, dup;
              beta = beamline->twiss0->betax;
              alpha = beamline->twiss0->alphax;
              du = (startingCoord[0] - trackingBuffer[0][0]) / beta;
              dup = (startingCoord[1] - trackingBuffer[0][1]) + alpha * du;
              dF = fabs(du) + fabs(dup);

              beta = beamline->twiss0->betay;
              alpha = beamline->twiss0->alphay;
              du = (startingCoord[2] - trackingBuffer[0][2]) / beta;
              dup = (startingCoord[3] - trackingBuffer[0][3]) + alpha * du;
              dF += fabs(du) + fabs(dup);
            }
            /* Log the data */
            if (!SDDS_SetRowValues(&SDDS_cmap, SDDS_SET_BY_INDEX | SDDS_PASS_BY_VALUE, ip,
                                   IC_X, x, IC_Y, y, IC_DELTA, delta,
                                   IC_SURVIVED, survived[0] * survived[1],
                                   IC_DF, dF, IC_LOGDF, log(dF + 1e-300),
                                   -1)) {
              SDDS_SetError("Problem setting SDDS row values (doChaosMap)");
              SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors | SDDS_EXIT_PrintErrors);
            }
          } else {
            /* Method using multi-turn tracking */
            /* Do the tracking */
            memcpy(trackingBuffer[0], startingCoord, sizeof(*startingCoord) * 6);
            trackingBuffer[0][6] = 1;
            p = run->p_central;
            memset(survived, 0, 3 * sizeof(*survived));
            if ((survived[0] = do_tracking(NULL, trackingBuffer, 1, NULL, beamline, &p, (double **)NULL, (BEAM_SUMS **)NULL, (long *)NULL,
                                           NULL, run, 0, TEST_PARTICLES, control->n_passes, 0,
                                           NULL, NULL, NULL, NULL, NULL))) {
              memcpy(trackingBuffer[1], startingCoord, sizeof(*startingCoord) * 6);
              trackingBuffer[1][0] += change_x;
              trackingBuffer[1][6] = 2;
              survived[1] = do_tracking(NULL, trackingBuffer + 1, 1, NULL, beamline, &p, (double **)NULL, (BEAM_SUMS **)NULL, (long *)NULL,
                                        NULL, run, 0, TEST_PARTICLES, control->n_passes, 0,
                                        NULL, NULL, NULL, NULL, NULL);

              memcpy(trackingBuffer[2], startingCoord, sizeof(*startingCoord) * 6);
              trackingBuffer[2][2] += change_y;
              trackingBuffer[2][6] = 3;
              survived[2] = do_tracking(NULL, trackingBuffer + 2, 1, NULL, beamline, &p, (double **)NULL, (BEAM_SUMS **)NULL, (long *)NULL,
                                        NULL, run, 0, TEST_PARTICLES, control->n_passes, 0,
                                        NULL, NULL, NULL, NULL, NULL);
            }
            /* Compute deltas */
            dJx = dJy = DBL_MAX;
            if (survived[0]) {
              if (survived[1]) {
                double J1, J2, beta, alpha, gamma;
                beta = beamline->twiss0->betax;
                alpha = beamline->twiss0->alphax;
                gamma = (1 + alpha * alpha) / beta;
                J1 = (sqr(trackingBuffer[0][0]) * gamma + 2 * alpha * trackingBuffer[0][0] * trackingBuffer[0][1] + sqr(trackingBuffer[0][1]) * beta) / 2;
                J2 = (sqr(trackingBuffer[1][0]) * gamma + 2 * alpha * trackingBuffer[1][0] * trackingBuffer[1][1] + sqr(trackingBuffer[1][1]) * beta) / 2;
                dJx = J2 - J1;
              }
              if (survived[2]) {
                double J1, J2, beta, alpha, gamma;
                beta = beamline->twiss0->betay;
                alpha = beamline->twiss0->alphay;
                gamma = (1 + alpha * alpha) / beta;
                J1 = (sqr(trackingBuffer[0][2]) * gamma + 2 * alpha * trackingBuffer[0][2] * trackingBuffer[0][3] + sqr(trackingBuffer[0][3]) * beta) / 2;
                J2 = (sqr(trackingBuffer[1][2]) * gamma + 2 * alpha * trackingBuffer[1][2] * trackingBuffer[1][3] + sqr(trackingBuffer[1][3]) * beta) / 2;
                dJy = J2 - J1;
              }
            }
            /* Log the data */
            if (!SDDS_SetRowValues(&SDDS_cmap, SDDS_SET_BY_INDEX | SDDS_PASS_BY_VALUE, ip,
                                   IC_X, x, IC_Y, y, IC_DELTA, delta,
                                   IC_DJX, dJx, IC_DJY, dJy, IC_SURVIVED, survived[0] * survived[1] * survived[2],
                                   IC_LOGDJX, log(fabs(dJx)), IC_LOGDJY, log(fabs(dJy)),
                                   -1)) {
              SDDS_SetError("Problem setting SDDS row values (doChaosMap)");
              SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors | SDDS_EXIT_PrintErrors);
            }
          }
          ip++;
          if (verbosity) {
#if USE_MPI
            if (myid == 1) {
              double newPercentage = 100 * (idelta * nx * ny + ix * ny + iy + 1.0) / (ndelta * nx * ny);
              if ((newPercentage - oldPercentage) >= 1) {
                printf("About %.1f%% done\n", newPercentage);
                oldPercentage = newPercentage;
                fflush(stdout);
              }
            }
#else
            printf("Done with particle %ld of %ld\n",
                   ix * ny * ndelta + iy * ndelta + idelta + 1, nx * ny * ndelta);
            fflush(stdout);
#endif
          }
        }
      }
    }
  }

  if (forward_backward > 0)
    free_beamlines(btBeamline);

  if (!inhibitFileSync)
    SDDS_DoFSync(&SDDS_cmap);
#if SDDS_MPI_IO
  if (!SDDS_MPI_WriteTable(&SDDS_cmap)) {
#else
  if (!SDDS_WriteTable(&SDDS_cmap)) {
#endif
    SDDS_SetError("Problem writing SDDS table (doChaosMap)");
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors | SDDS_EXIT_PrintErrors);
  }

#if USE_MPI
  /* disable output from first slave */
  if (myid == 1) {
#  if defined(_WIN32)
    freopen("NUL", "w", stdout);
#  else
    if (!freopen("/dev/null", "w", stdout)) {
      perror("freopen failed");
      exit(EXIT_FAILURE);
    }
#  endif
  }
#endif

#if USE_MPI
  MPI_Barrier(MPI_COMM_WORLD);
#endif
  return (1);
}

void finishChaosMap() {
  if (SDDS_IsActive(&SDDS_cmap) && !SDDS_Terminate(&SDDS_cmap)) {
    SDDS_SetError("Problem terminating SDDS output (finishChaosMap)");
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors | SDDS_EXIT_PrintErrors);
  }
}
