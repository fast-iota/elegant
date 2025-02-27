#include "mdb.h"
#include "track.h"

char **addPatterns(long *patterns, char *input0) {
  char *input, *ptr, **pattern;
  *patterns = 0;
  if (input0 == NULL)
    return NULL;
  cp_str(&input, input0);
  pattern = NULL;
  while ((ptr = get_token(input))) {
    pattern = SDDS_Realloc(pattern, sizeof(*pattern) * (*patterns + 1));
    cp_str(&(pattern[*patterns]), ptr);
    if (has_wildcards(pattern[*patterns]) && strchr(pattern[*patterns], '-'))
      pattern[*patterns] = expand_ranges(pattern[*patterns]);
    *patterns += 1;
  }
  free(input);
  return pattern;
}

int matchesPatternList(char **pattern, long patterns, char *input) {
  long i;
  for (i = 0; i < patterns; i++) {
    if (wild_match(input, pattern[i])) {
      return 1;
    }
  }
  return 0;
}
