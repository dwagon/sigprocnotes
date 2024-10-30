/*
 * texproc: embedded command preprocessor for TeX and LaTeX
 * (c) 1992-2024 Dougal Scott
 * Any comments, criticisms, patches to 
 * dougal.scott@gmail.com
 *
 * Convert LaTeX:
 * ....
 * %# gnuplot
 * plot sin, cos, and tan
 * %#
 * ....
 *
 * to 
 *
 * ....
 * \begin{picture}
 * \lotsadots
 * \end{picture}
 * ....
 */

#define TRACE(x) /* x */

#include <stdio.h>
#include <stdlib.h>
#include <strings.h>
#include <sys/types.h>
#include <unistd.h>
#include <errno.h>
#include <limits.h>

extern int errno;
FILE * tmp;

/****************************************************************************************************************/
void process(char * buff) {
  char tmpname[BUFSIZ], /* Name of tmp file */
    progname[BUFSIZ], /* Name of program to execute */
    cmdline[BUFSIZ], /* What to feed to popen */
    outbuff[BUFSIZ]; /* What to take output of program from */
  FILE * p;

  strncpy(tmpname, "/tmp/texproc.out", BUFSIZ);
  if ((tmp = fopen(tmpname, "w")) == NULL) {
    fprintf(stderr, "Could not open %s for writing\n", tmpname);
    fprintf(stderr, "Program aborting\n");
    exit(-1);
  }
  TRACE(fprintf(stderr, "Saving to tmp file %s\n", tmpname));
  strncpy(progname, & buff[3], BUFSIZ);  // Start from after the initial #%
  progname[strcspn(progname, "\n")] = '\0';
  fprintf(stdout, "%% Including output from %s\n", progname);
  /* Put buffer to file for executing */
  while (fgets(buff, BUFSIZ, stdin) != NULL) {
    if (buff[0] == '%' && buff[1] == '#') {
      fclose(tmp);
      buff[strcspn(buff, "\n")] = '\0';
      TRACE(fprintf(stderr, "Executing >%s< >%s<\n", progname, tmpname));
      sprintf(cmdline, "%s %s", progname, tmpname);
      fprintf(stderr, "%s\n", progname);
      if ((p = popen(cmdline, "r")) == NULL) {
        fprintf(stderr, "Could not open pipe to %s\n", cmdline);
        exit(-1);
      }
      while (fgets(outbuff, BUFSIZ, p) != NULL)
        fprintf(stdout, "%s", outbuff);
      pclose(p);
      unlink(tmpname);
      return;
    } else {
      TRACE(fprintf(stderr, "%s", buff));
      fprintf(tmp, "%s", buff);
    }
  }
}

/****************************************************************************************************************/
int main(int argc, char ** argv) {
  char buff[BUFSIZ];

  while (fgets(buff, BUFSIZ, stdin) != NULL) {
    if (buff[0] == '%' && buff[1] == '#')
      process(buff);
    else
      fprintf(stdout, "%s", buff);
  }
  return (0);
}