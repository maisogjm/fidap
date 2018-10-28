#include <stdio.h>
#include <stdlib.h>
#include <math.h>
float F2(float F_test);
float betai(float a, float b, float x);
float golden(float ax, float bx, float cx, float (*f)(float), float tol,
	float *xmin);

float v1, v2;
float prob_thresh;
 
main(int argc, char *argv[])
{
    float F;
    float min;

/* ---------------------------------------------------------------------------- */
/*  CHECK TO SEE IF CORRECT NUMBER OF ARGUMENTS WERE PASSED.                    */

    if(argc<4) {
	printf("F: need probability threshold and v1 and v2.\n");
	exit(1);
    }

/* ---------------------------------------------------------------------------- */
/*  READ IN COMMAND LINE ARGUMENTS.                                             */

    prob_thresh = atof(argv[1]);
    if ((prob_thresh>1.) || (prob_thresh<0.)) {
        printf("Probability threshold must be between 0.0 and 1.0.\n");
        exit(2);
    }
    v1     = atof(argv[2]);
    v2     = atof(argv[3]);
/*    printf("prob_thresh = %f\n",prob_thresh);
    printf("df          = %4.0f\n",df); */

    min=golden(0,3,50,(float (*) (float)) F2,0.0001,&F);

/*    printf("For p=%f, df=%d, threshold 2-tailed t=%f\n",prob_thresh,(int)df,t);
    printf("Deviation = %f\n",min); */
    printf("%g\n",F);
}

/* ---------------------------------------------------------------------------- */
/* F distribution.								*/
float F2(float F_test)
{
/*    printf("prob_thresh = %f, current t_prob = %f, t = %f\n",prob_thresh,
	1.-betai(df/2.,0.5,(df/(df+(F_test*F_test)))),F_test);
    printf("%f\n",(float) fabs((double) (prob_thresh-betai(df/2.,0.5,(df/(df+(F_test*F_test))))))); */
    return (float) fabs(prob_thresh-betai(v2/2.,v1/2,(v2/(v2+(v1*F_test)))));
}
