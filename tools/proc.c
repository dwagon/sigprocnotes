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

void process(char * buff) {
  char tmpname[80], /* Name of tmp file */
    progname[80], /* Name of program to execute */
    cmdline[80], /* What to feed to popen */
    outbuff[BUFSIZ]; /* What to take output of prgram from */
  FILE * p;

  strcpy(tmpname, "/tmp/PrcXXXXXX");
  mkstemp(tmpname);
  if ((tmp = fopen(tmpname, "w")) == NULL) {
    fprintf(stderr, "Could not open %s for writing\n", tmpname);
    fprintf(stderr, "Program aborting\n");
    exit(-1);
  }
  TRACE(fprintf(stderr, "Saving to tmp file %s\n", tmpname));
  strcpy(progname, & buff[3]);
  fprintf(stdout, "%% Including output from %s\n", progname);
  /* Put buffer to file for executing */
  while (fgets(buff, BUFSIZ, stdin) != NULL) {
    if (buff[0] == '%' && buff[1] == '#') {
      fclose(tmp);
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
      TRACE(fprintf(stderr, "%s\n", buff));
      fprintf(tmp, "%s\n", buff);
    }
  }
}

int main(int argc, char ** argv) {
  char buff[BUFSIZ];

  while (fgets(buff, BUFSIZ, stdin) != NULL) {
    if (buff[0] == '%' && buff[1] == '#')
      process(buff);
    else
      fprintf(stdout, "%s\n", buff);
  }
  return (0);
}