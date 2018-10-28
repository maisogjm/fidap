#include <stdio.h>
#include <stdlib.h>
#include <math.h>
main(int argc, char *argv[])
{
    float a;

    /* **************************************************************** */
    /* VERIFY NUMBER OF COMMAND LINE ARGUMENTS.				*/
    if (argc < 2) {
	printf("sqrt: need a non-negative number.\n");
	printf("      Will return square root of that number.\n");
	exit(1);
    }

    /* **************************************************************** */
    /* READ IN NUMBER.							*/
    a = atof(argv[1]);

    /* **************************************************************** */
    /* MAKE SURE a IS NON-NEGATIVE.					*/
    if (a<0) {
	printf("sqrt: number must be non-negative.\n");
	exit(2);
    }

    /* **************************************************************** */
    /* CALCULATE SQUARE ROOT.                                           */
    printf("%g\n",(float)sqrt((double)a));
}
