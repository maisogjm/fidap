#include <stdio.h>
#include <stdlib.h>

main (int argc, char *argv[])
{
    /* 09.08.93 : JM Maisog, M.D.
	Skips a specified number of bytes from an input file, and
	then passes a specified number of bytes to an output file,
	APPENDING to that file if it already exist.		   */

    FILE *fIN, *fOUT;
    int c;
    int count;

    if (argc < 5) {
	printf("excise: need input file, output file, number of\n");
	printf("        bytes to skip, number of bytes to keep.\n");
	exit(1);
    }

    /* OPEN FILES FOR INPUT AND OUTPUT.					*/
    if ((fIN=fopen(argv[1],"r")) == NULL) {
	printf("excise: can't open %s\n",argv[1]);
	exit(2);
    }
    if ((fOUT=fopen(argv[2],"a")) == NULL) {
	printf("excise: can't open %s\n",argv[2]);
	exit(3);
    }

    /* SKIP NUMBER OF BYTES SPECIFIED.					*/
    count = 0;
    while((count<atoi(argv[3])) && ((c=fgetc(fIN)) != EOF))
	count++;

    /* SPIT OUT NUMBER OF BYTES SPECIFIED.				*/
    count = 0;
    while((c=fgetc(fIN)) != EOF && (count<atoi(argv[4]))) {
	fputc(c,fOUT);
	count++;
    }

    fclose(fIN);
    fclose(fOUT);
}
