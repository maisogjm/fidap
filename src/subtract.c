#include <stdio.h>
#include <stdlib.h>
main(int argc, char *argv[])
{
    float a,b;

    /* **************************************************************** */
    /* VERIFY NUMBER OF COMMAND LINE ARGUMENTS.				*/
    if (argc < 3) {
	printf("subtract: need two numbers.\n");
	printf("          Will return the first number minus the second.\n");
	exit(1);
    }

    /* **************************************************************** */
    /* READ IN TWO NUMBERS.						*/
    a = atof(argv[1]);
    b = atof(argv[2]);

    /* **************************************************************** */
    /* DO DIVISION.							*/
    printf("%g\n",a-b);
}
