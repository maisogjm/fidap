#include <stdio.h>
#include <stdlib.h>
main(int argc, char *argv[])
{
    float a,b;

    /* **************************************************************** */
    /* VERIFY NUMBER OF COMMAND LINE ARGUMENTS.				*/
    if (argc < 3) {
	printf("div: need two numbers.\n");
	printf("     Will return the first number divided by the second.\n");
	exit(1);
    }

    /* **************************************************************** */
    /* READ IN TWO NUMBERS.						*/
    a = atof(argv[1]);
    b = atof(argv[2]);

    /* **************************************************************** */
    /* MAKE SURE b IS NOT ZERO.						*/
    if (b==0) {
	printf("div: divide by zero error.\n");
	exit(2);
    }

    /* **************************************************************** */
    /* DO DIVISION.							*/
    printf("%g\n",a/b);
    exit(0);
}
