#include <stdio.h>
#include <stdlib.h>
#include <math.h>
float r2(float r_test);
float betai(float a, float b, float x);
float golden(float ax, float bx, float cx, float (*f)(float), float tol,
	float *xmin);

float df;
float prob_thresh;
 
main(int argc, char *argv[])
{
    float r;
    float min;

/* ---------------------------------------------------------------------------- */
/*  CHECK TO SEE IF CORRECT NUMBER OF ARGUMENTS WERE PASSED.                    */

    if(argc<3) {
	printf("r: need probability threshold and df (usually N-2).\n");
	exit(1);
    }

/* ---------------------------------------------------------------------------- */
/*  READ IN COMMAND LINE ARGUMENTS.                                             */

    prob_thresh = atof(argv[1]);
    df          = atof(argv[2]);

    min=golden(0.,0.5,1.,(float (*) (float)) r2,0.0001,&r);

    printf("%g\n",r);
}

/* ---------------------------------------------------------------------------- */
/* r distribution.								*/
float r2(float r_test)
{
    float t_test;
    t_test=(float)(sqrt((double)(df))*r_test)/(float)(sqrt((double)(1.-r_test*r_test)));
    return (float) fabs(prob_thresh-betai(df/2.,0.5,(df/(df+(t_test*t_test)))));
}
