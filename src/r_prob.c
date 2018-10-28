#include <stdio.h>
#include <stdlib.h>
#include <math.h>
float r2(float r);
float betai(float a, float b, float x);
float golden(float ax, float bx, float cx, float (*f)(float), float tol,
	float *xmin);

float df;
float r;
 
main(int argc, char *argv[])
{
    float r,t_test;
    float min;

/* ---------------------------------------------------------------------------- */
/*  CHECK TO SEE IF CORRECT NUMBER OF ARGUMENTS WERE PASSED.                    */

    if(argc<3) {
	printf("r_prob: need test statistic and df (usually N-2).\n");
	exit(1);
    }

/* ---------------------------------------------------------------------------- */
/*  READ IN COMMAND LINE ARGUMENTS.                                             */

    r  = atof(argv[1]);
    df = atof(argv[2]);

    t_test=(float)(sqrt((double)(df))*r)/(float)(sqrt((double)(1.-r*r)));
    printf("%g\n",betai(df/2.,0.5,(df/(df+(t_test*t_test)))));
}
