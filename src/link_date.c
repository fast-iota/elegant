/*************************************************************************\
* Copyright (c) 2002 The University of Chicago, as Operator of Argonne
* National Laboratory.
* Copyright (c) 2002 The Regents of the University of California, as
* Operator of Los Alamos National Laboratory.
* This file is distributed subject to a Software License Agreement found
* in the file LICENSE that is included with this distribution. 
\*************************************************************************/

#include <stdio.h>
#if USE_MPI
#  include "mpi.h"
#endif

void link_date(void) {
#if USE_MPI
#  ifdef MPI_MAX_LIBRARY_VERSION_STRING
  int len;
  char mpi_lib_ver[MPI_MAX_LIBRARY_VERSION_STRING];
#  endif
#endif
#ifdef SVN_VERSION
  printf("Link date: %s %s, SVN revision: %s\n", __DATE__, __TIME__, SVN_VERSION);
#else
  printf("Link date: %s %s\n", __DATE__, __TIME__);
#endif

#if USE_MPI
#  ifdef MPI_MAX_LIBRARY_VERSION_STRING
  MPI_Get_library_version(mpi_lib_ver, &len);
  printf("MPI library version: %s\n", mpi_lib_ver);
#  endif
#endif
  fflush(stdout);
}
