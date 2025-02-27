/*************************************************************************\
 * Copyright (c) 2002 The University of Chicago, as Operator of Argonne
 * National Laboratory.
 * Copyright (c) 2002 The Regents of the University of California, as
 * Operator of Los Alamos National Laboratory.
 * This file is distributed subject to a Software License Agreement found
 * in the file LICENSE that is included with this distribution. 
\*************************************************************************/

#include "mdb.h"
//#include "mdbsun.h"
#include "track.h"


// TODO: is this a bug? flags not set
TRACKING_CONTEXT trackingContext =
  {"", -1, 0, -1, NULL, NULL, 0.0, 0.0, "", 0
#if USE_MPI
   ,
   -1
#endif
};

void getTrackingContext(TRACKING_CONTEXT *trackingContext0) {
  memcpy(trackingContext0, &trackingContext, sizeof(trackingContext));
}
#if TURBO_STRINGS
static char nullstr[] = "\0";
#endif
void setTrackingContext(char *name, long occurence, long type, char *rootname, ELEMENT_LIST *eptr) {
#if USE_MPI
  trackingContext.myid = myid;
#endif
  trackingContext.element = eptr;
  trackingContext.sliceAnalysis = NULL;
  trackingContext.zStart = 0;
  trackingContext.zEnd = 0;
  trackingContext.step = 0;
  trackingContext.flags = 0;

#if TURBO_STRINGS
  if (name)
    trackingContext.elementName = name;
  else
    trackingContext.elementName = nullstr;
#else
  if (name)
    strncpy(trackingContext.elementName, name, CONTEXT_BUFSIZE);
  else
    trackingContext.elementName[0] = 0;
#endif

  trackingContext.elementOccurrence = occurence;
  trackingContext.elementType = type;
#if TURBO_STRINGS
  if (rootname)
    trackingContext.rootname = rootname;
  else
    trackingContext.rootname = nullstr;
#else
  if (rootname)
    strncpy(trackingContext.rootname, rootname, CONTEXT_BUFSIZE);
  else
    trackingContext.rootname[0] = 0;
#endif

}
