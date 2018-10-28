#include <math.h>
#include <stdio.h>
#include <stdlib.h>
float golden(float ax, float bx, float cx, float (*f)(float), float tol,
        float *xmin);
float gammq(float a, float x);
float chi2(float chi2_test);

    float prob_thresh,        /* Probability threshold.                         */
              df;        /* degrees of freedom.                                 */
main(int argc, char *argv[])
{
    float min,chi2_test;

/* ---------------------------------------------------------------------------- */
/*  CHECK TO SEE IF CORRECT NUMBER OF ARGUMENTS WERE PASSED.                    */

    if(argc<3) {
        printf("chi2: need probability and degrees of freedom.\n");
        exit(1);
    }

/* ---------------------------------------------------------------------------- */
/*  READ IN COMMAND LINE ARGUMENTS.                                             */

    prob_thresh = atof(argv[1]);
    df     = atof(argv[2]);

/* ---------------------------------------------------------------------------- */

    min=golden(0.,5,10000,(float (*) (float)) chi2,0.0001,&chi2_test);
    printf("%g\n",chi2_test);
}
/* ---------------------------------------------------------------------------- */
/* chi2 distribution.                                                           */
float chi2(float chi2_test)
{
    return (float) fabs((double)(prob_thresh-gammq(df/2.,chi2_test/2.)));
}

