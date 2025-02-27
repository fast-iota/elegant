/* This is not used by elegant.  It is used by other programs that use some of elegant's
 * subroutines.  Elegant.c has a copy that is used by elegant itself.
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>

void bombElegant(char *error, char *usage) {
  if (error)
    fprintf(stderr, "error: %s\n", error);
  if (usage)
    fprintf(stderr, "usage: %s\n", usage);
  exit(1);
}

void exitElegant(long status) {
  exit(status);
}

void bombElegantVA(char *template, ...) {
  char *p;
  char c, *s;
  int i;
  long j;
  va_list argp;
  double d;

  va_start(argp, template);
  p = template;
  while (*p) {
    if (*p == '%') {
      switch (*++p) {
      case 'l':
        switch (*++p) {
        case 'd':
          j = va_arg(argp, long int);
          printf("%ld", j);
          break;
        case 'e':
          d = va_arg(argp, double);
          printf("%21.15le", d);
          break;
        case 'f':
          d = va_arg(argp, double);
          printf("%lf", d);
          break;
        case 'g':
          d = va_arg(argp, double);
          printf("%21.15lg", d);
          break;
        default:
          printf("%%l%c", *p);
          break;
        }
        break;
      case 'c':
        c = va_arg(argp, int);
        putchar(c);
        break;
      case 'd':
        i = va_arg(argp, int);
        printf("%d", i);
        break;
      case 's':
        s = va_arg(argp, char *);
        fputs(s, stdout);
        break;
      case 'e':
        d = va_arg(argp, double);
        printf("%21.15e", d);
        break;
      case 'f':
        d = va_arg(argp, double);
        printf("%f", d);
        break;
      case 'g':
        d = va_arg(argp, double);
        printf("%21.15g", d);
        break;
      default:
        printf("%%%c", *p);
        break;
      }
    } else {
      putchar(*p);
    }
    p++;
  }
  exit(1);
}
