#include <stdio.h>
#include <stdlib.h>
float betai(float a, float b, float x);

main(int argc, char *argv[])
{
    float t_test,	/* t-test						*/
	      df;	/* degrees of freedom.					*/

/* ---------------------------------------------------------------------------- */
/*  CHECK TO SEE IF CORRECT NUMBER OF ARGUMENTS WERE PASSED.                    */

    if(argc<3) {
	printf("t_prob: need test statistic and degrees of freedom.\n");
        printf("        Will return TWO-sided significance.\n");
	exit(1);
    }

/* ---------------------------------------------------------------------------- */
/*  READ IN COMMAND LINE ARGUMENTS.						*/

    t_test = atof(argv[1]);
    df     = atof(argv[2]);
/*    printf("t_test = %f\n",t_test);
    printf("df     = %4.0f\n",df); */

/* ---------------------------------------------------------------------------- */

    printf("%g\n",betai(df/2.,0.5,(df/(df+(t_test*t_test)))));
}
