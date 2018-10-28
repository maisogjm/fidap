#include <stdio.h>
#include <stdlib.h>
#include <math.h>
main(int argc, char *argv[])
{
    float a;

    /* **************************************************************** */
    /* VERIFY NUMBER OF COMMAND LINE ARGUMENTS.				*/
    if (argc < 2) {
	printf("loge: need a positive number.\n");
	printf("      Will return natural log of that number.\n");
	exit(1);
    }

    /* **************************************************************** */
    /* READ IN NUMBER.							*/
    a = atof(argv[1]);

    /* **************************************************************** */
    /* MAKE SURE a IS POSITIVE.     				*/
    if (a<=0) {
	printf("loge: number must be positive.\n");
	exit(2);
    }

    /* **************************************************************** */
    /* CALCULATE NATURAL LOG.                                           */
    printf("%g\n",(float)log((double)a));
}
