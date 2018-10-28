#include <stdio.h>
#include <stdlib.h>
#include <math.h>
main(int argc, char *argv[])
{
    float a;

    /* **************************************************************** */
    /* VERIFY NUMBER OF COMMAND LINE ARGUMENTS.				*/
    if (argc < 2) {
	printf("log10: need a positive number.\n");
	printf("       Will return log base 10 of that number.\n");
	exit(1);
    }

    /* **************************************************************** */
    /* READ IN NUMBER.							*/
    a = atof(argv[1]);

    /* **************************************************************** */
    /* MAKE SURE a IS POSITIVE.     				*/
    if (a<=0) {
	printf("log10: number must be positive.\n");
	exit(2);
    }

    /* **************************************************************** */
    /* CALCULATE LOG.                                                   */
    printf("%g\n",(float)log((double)a)/(float)(log(10.)));
}
