#include <stdio.h>
#include <stdlib.h>
#include <math.h>
main(int argc, char *argv[])
{
    float a,cum_prod;
    int index;

    /* **************************************************************** */
    /* VERIFY NUMBER OF COMMAND LINE ARGUMENTS.				*/
    if (argc < 3) {
	printf("geomean: need two or more numbers.\n");
	exit(1);
    }

    /* **************************************************************** */
    /* LOOP OVER COMMAND LINE ARGUMENTS AND COMPUTE GEOMETRIC MEAN.     */

    cum_prod=1.;
    for (index=1; index<argc; index++)
        cum_prod*=atof(argv[index]);

    /* **************************************************************** */
    /* DO DIVISION.							*/
    printf("%g\n",pow(cum_prod,1./((float)argc-1.)));
}
