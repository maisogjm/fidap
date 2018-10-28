#include <stdio.h>
#include <stdlib.h>
#include <math.h>

main(int argc, char *argv[])
{
    float a,b;          /* Lower and upper bounds of range over which normal    */
                        /* curve will be integrated.                            */
    float SD;           /* Standard deviation of gaussian curve to be integrated*/
    float ca,ya,za;
    float area_a;


/* ---------------------------------------------------------------------------- */
/*  CHECK TO SEE IF CORRECT NUMBER OF ARGUMENTS WERE PASSED.                    */

    if (argc!=2) {
        printf("z_prob: need Z-score.\n");
        printf("Will return area under normal curve from Z-score to +Inf.\n");
        exit(1);
    }

/* ---------------------------------------------------------------------------- */
/* READ IN a.                                                                   */

    a  = atof(argv[1]);

/* ---------------------------------------------------------------------------- */
/* CALCULATE y, z & c.                                                          */

    ya=1./(1.+0.2316419*fabs(a));
    za=0.3989423*exp(-(a*a)/2.);
    ca=1-za*(1.330274*ya*ya*ya*ya*ya-1.821256*ya*ya*ya*ya+1.781478*ya*ya*ya-0.356538*ya*ya+.3193815*ya);
    if (a>=0)
        area_a=1.-ca;
    else
        area_a=ca;

    printf("%g\n",area_a);
}

