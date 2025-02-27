#if defined(_WIN32)
#  ifdef __cplusplus
extern "C" {
#  endif
short has_occurence_string(char *string);
#  ifdef __cplusplus
}
#  endif

#  include <regex>

#  if defined(_MINGW)
short has_occurence_string(char *string) {
  const std::regex pattern(".*%[0-9]*ld.*");
  return std::regex_match(string, pattern);
}
#  else
short has_occurence_string(char *string) {
  const std::tr1::regex pattern(".*%[0-9]*ld.*");
  return std::tr1::regex_match(string, pattern);
}
#  endif
#endif
