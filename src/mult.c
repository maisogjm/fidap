#include <stdio.h>
#include <stdlib.h>
main(int argc, char *argv[])
{
    float a,cum_prod;
    int index;

    /* **************************************************************** */
    /* VERIFY NUMBER OF COMMAND LINE ARGUMENTS.                         */
    if (argc < 3) {
        printf("mult: need two or more numbers.\n");
        printf("      Will return the product of those numbers.\n");
        exit(1);
    }

    /* **************************************************************** */
    /* LOOP OVER COMMAND LINE ARGUMENTS AND COMPUTE PRODUCT.            */

    cum_prod=1.;
    for (index=1; index<argc; index++)
        cum_prod*=atof(argv[index]);

    /* **************************************************************** */
    /* PRINTOUT RESULT.                                                 */
    printf("%g\n",cum_prod);
}
