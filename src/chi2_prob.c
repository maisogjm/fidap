#include <stdio.h>
#include <stdlib.h>
float gammq(float a, float x);

main(int argc, char *argv[])
{
    float chi2_test,	/* t-test						*/
	      df;	/* degrees of freedom.					*/

/* ---------------------------------------------------------------------------- */
/*  CHECK TO SEE IF CORRECT NUMBER OF ARGUMENTS WERE PASSED.                    */

    if(argc<3) {
	printf("chi2_prob: need test statistic and degrees of freedom.\n");
	exit(1);
    }

/* ---------------------------------------------------------------------------- */
/*  READ IN COMMAND LINE ARGUMENTS.						*/

    chi2_test = atof(argv[1]);
    df     = atof(argv[2]);
/*    printf("chi2_test = %f\n",chi2_test);
    printf("df     = %4.0f\n",df); */

/* ---------------------------------------------------------------------------- */

    printf("%g\n",gammq(df/2.,chi2_test/2.));
}
