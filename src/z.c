#include <stdio.h>
#include <stdlib.h>
#include <math.h>
float golden(float ax, float bx, float cx, float (*f)(float), float tol,
        float *xmin);
float prob_thresh;
float normal(float zscore);

main(int argc, char *argv[])
{
    float min,z;

/* ---------------------------------------------------------------------------- */
/*  CHECK TO SEE IF CORRECT NUMBER OF ARGUMENTS WERE PASSED.                    */

    if (argc!=2) {
        printf("z: need probability.\n");
        printf("Will return Z-score with given 1-tailed significance.\n");
        exit(1);
    }

/* ---------------------------------------------------------------------------- */
/* READ IN prob_thresh.                                                         */

    prob_thresh  = atof(argv[1]);
    if ((prob_thresh < 0.) || (prob_thresh > 1.)) {
	fprintf(stderr,"Probability must be between 0 and 1!\n");
	exit(2);
    }

/* ---------------------------------------------------------------------------- */
/* CALCULATE y, z & c.                                                          */

    min=golden(-10.,0.,10.,(float (*) (float)) normal,0.0001,&z);
    printf("%g\n",z);
    exit(0);
}
/* ---------------------------------------------------------------------------- */
/* Normal distribution.                                                         */
float normal(float zscore)
{
    float area_a,ya,za,ca;

    ya=1./(1.+0.2316419*fabs(zscore));
    za=0.3989423*exp(-(zscore*zscore)/2.);
    ca=1-za*(1.330274*ya*ya*ya*ya*ya-1.821256*ya*ya*ya*ya+1.781478*ya*ya*ya-0.356538*ya*ya+.3193815*ya);
    if (zscore>=0)
        area_a=1.-ca;
    else
        area_a=ca;
    return (float) fabs(prob_thresh-area_a);
}
