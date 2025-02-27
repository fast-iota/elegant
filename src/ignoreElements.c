/*************************************************************************\
* Copyright (c) 2017 The University of Chicago, as Operator of Argonne
* National Laboratory.
* This file is distributed subject to a Software License Agreement found
* in the file LICENSE that is included with this distribution. 
\*************************************************************************/

#include "mdb.h"
#include "track.h"
#include "match_string.h"

#define DEBUG 0

typedef struct {
  char *name, *type, *exclude;
  long completely;
} IGNORE_ELEMENT_SPEC;

static IGNORE_ELEMENT_SPEC *ignoreElementsSpec = NULL;
static long ignoreElementsSpecs = 0;

long countIgnoreElementsSpecs(long completely0) {
  long count = 0;
  long i;
  for (i = 0; i < ignoreElementsSpecs; i++) {
    if (ignoreElementsSpec[i].completely == completely0)
      count++;
  }
  return count;
}

void addIgnoreElementsSpec(char *name, char *type, char *exclude, long completely0) {
  if (!(ignoreElementsSpec = SDDS_Realloc(ignoreElementsSpec,
                                          sizeof(*ignoreElementsSpec) * (ignoreElementsSpecs + 1))))
    bombElegant("memory allocation failure", NULL);
  ignoreElementsSpec[ignoreElementsSpecs].name = NULL;
  ignoreElementsSpec[ignoreElementsSpecs].type = NULL;
  ignoreElementsSpec[ignoreElementsSpecs].exclude = NULL;
  ignoreElementsSpec[ignoreElementsSpecs].completely = completely0;
  if ((name &&
       !SDDS_CopyString(&ignoreElementsSpec[ignoreElementsSpecs].name, name)) ||
      (type &&
       !SDDS_CopyString(&ignoreElementsSpec[ignoreElementsSpecs].type, type)) ||
      (exclude &&
       !SDDS_CopyString(&ignoreElementsSpec[ignoreElementsSpecs].exclude, exclude)))
    bombElegant("memory allocation failure", NULL);

  ignoreElementsSpecs++;
}

void clearIgnoreElementsSpecs() {
  while (ignoreElementsSpecs--) {
    if (ignoreElementsSpec[ignoreElementsSpecs].name)
      free(ignoreElementsSpec[ignoreElementsSpecs].name);
    if (ignoreElementsSpec[ignoreElementsSpecs].type)
      free(ignoreElementsSpec[ignoreElementsSpecs].type);
    if (ignoreElementsSpec[ignoreElementsSpecs].exclude)
      free(ignoreElementsSpec[ignoreElementsSpecs].exclude);
  }
  free(ignoreElementsSpec);
  ignoreElementsSpec = NULL;
}

long ignoreElement(char *name, long type, long completelyOnly) {
  long i;
  for (i = 0; i < ignoreElementsSpecs; i++) {
    if (ignoreElementsSpec[i].exclude && wild_match(name, ignoreElementsSpec[i].exclude))
      continue;
    if (ignoreElementsSpec[i].name && !wild_match(name, ignoreElementsSpec[i].name))
      continue;
    if (ignoreElementsSpec[i].type && !wild_match(entity_name[type], ignoreElementsSpec[i].type))
      continue;
    if (ignoreElementsSpec[i].completely != completelyOnly)
      continue;
    return 1;
  }
  return 0;
}

#include "ignoreElements.h"

void setupIgnoreElements(NAMELIST_TEXT *nltext, RUN *run,
                         LINE_LIST *beamline) {
  long i;
  /* process the namelist text */
  set_namelist_processing_flags(STICKY_NAMELIST_DEFAULTS);
  set_print_namelist_flags(0);
  if (processNamelist(&ignore_elements, nltext) == NAMELIST_ERROR)
    bombElegant(NULL, NULL);
  if (echoNamelists)
    print_namelist(stdout, &ignore_elements);

  if (clear_all) {
    clearIgnoreElementsSpecs();
    if (!name && !type)
      return;
  }
  if (disable)
    return;

  if (!name || !strlen(name))
    bombElegant("no name given", NULL);
  str_toupper(name);
  if (has_wildcards(name) && strchr(name, '-'))
    name = expand_ranges(name);
  if (type) {
    str_toupper(type);
    if (has_wildcards(type) && strchr(type, '-'))
      type = expand_ranges(type);
    for (i = 0; i < N_TYPES; i++)
      if (wild_match(entity_name[i], type))
        break;
    if (i == N_TYPES) {
      fprintf(stderr, "type pattern %s does not match any known type", type);
      exitElegant(1);
    }
  }
  if (exclude) {
    str_toupper(exclude);
    if (has_wildcards(exclude) && strchr(exclude, '-'))
      exclude = expand_ranges(exclude);
  }

  addIgnoreElementsSpec(name, type, exclude, completely);
}
