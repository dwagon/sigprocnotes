/*
 * Convert LaTeX:
 * blah, blah, blah, ....
 * %# gnuplot
 * plot sin, cos, and tan
 * %#
 * blah, blah, blah, ....
 *
 * to 
 * blah, blah, blah, ....
 * \begin{picture}
 * \lotsadots
 * \end{picture}
 * blah, blah, blah, ....
 */

#define TRACE(x)	/* x */

#include <stdio.h>
#include <strings.h>
#include <sys/types.h>
#include <unistd.h>
#include <errno.h>

extern int errno;
FILE *tmp;

main(argc,argv)
int argc;
char **argv;
{
char buff[BUFSIZ];

while(gets(buff)!=NULL) {
	if(buff[0]=='%' && buff[1]=='#')
		process(buff);
	else
		fprintf(stdout,"%s\n",buff);
	}
return(0);
}

process(buff)
char *buff;
{
char 	tmpname[80],		/* Name of tmp file */
		progname[80];		/* Name of program to execute */	
pid_t pid;
int ret;

strcpy(tmpname,"/tmp/PrcXXXXXX");
mktemp(tmpname);
if((tmp=fopen(tmpname,"w"))==NULL) {
	fprintf(stderr,"Could not open %s for writing\n",tmpname);
	fprintf(stderr,"Program aborting\n");
	exit(-1);
	}
TRACE(fprintf(stderr,"Saving to tmp file %s\n",tmpname));
strcpy(progname,&buff[3]);
fprintf(stdout,"%% Including output from %s\n",progname);
/* Put buffer to file for executing */
while(gets(buff)!=NULL) {
	if(buff[0]=='%' && buff[1]=='#') {
		fclose(tmp);
		TRACE(fprintf(stderr,"forking\n"));
		pid=fork();
		TRACE(fprintf(stderr,"PID=%d\n",pid));
		if(pid==(pid_t)0) {
			TRACE(fprintf(stderr,"Forked:%s,%s\n",progname,tmpname));
			ret=execlp(progname,tmpname,NULL);
			fprintf(stderr,"Could not execute %s. Error=%d\n",progname,ret);
			fprintf(stderr,"Program aborting\n");
			exit(-1);
			}
		unlink(tmpname);
		return(0);
		}
	else {
		TRACE(fprintf(stderr,"%s\n",buff));
		fprintf(tmp,"%s\n",buff);
		}
	}
}
