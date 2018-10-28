#include <math.h>
#include <stdio.h>
#include <stdlib.h>
float betai(float a, float b, float x);

main(int argc, char *argv[])
{
    float F_test,	/* F-test						*/
	      v1,v2;	/* degrees of freedom.					*/

/* ---------------------------------------------------------------------------- */
/*  CHECK TO SEE IF CORRECT NUMBER OF ARGUMENTS WERE PASSED.                    */

    if(argc<4) {
	printf("F_prob: need test statistic and v1 and v2.\n");
	exit(1);
    }

/* ---------------------------------------------------------------------------- */
/*  READ IN COMMAND LINE ARGUMENTS.						*/

    F_test = (float)fabs((double)atof(argv[1]));
    v1     = atof(argv[2]);
    v2     = atof(argv[3]);
/*    printf("F_test = %f\n",F_test);
    printf("df     = %4.0f\n",df); */

/* ---------------------------------------------------------------------------- */

    printf("%g\n",betai(v2/2.,v1/2,(v2/(v2+(v1*F_test)))));
}
