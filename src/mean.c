#include <stdio.h>
#include <stdlib.h>
main(int argc, char *argv[])
{
    float a,cum_sum;
    int index;

    /* **************************************************************** */
    /* VERIFY NUMBER OF COMMAND LINE ARGUMENTS.				*/
    if (argc < 3) {
	printf("mean: need two or more numbers.\n");
	exit(1);
    }

    /* **************************************************************** */
    /* LOOP OVER COMMAND LINE ARGUMENTS AND COMPUTE PRODUCT.            */

    cum_sum=0.;
    for (index=1; index<argc; index++)
        cum_sum+=atof(argv[index]);

    /* **************************************************************** */
    /* CALCULATE AVERAGE.	*/
    printf("%g\n",cum_sum/(float)(argc-1));
}
