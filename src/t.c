#include <stdio.h>
#include <stdlib.h>
#include <math.h>
float t2(float t_test);
float betai(float a, float b, float x);
float golden(float ax, float bx, float cx, float (*f)(float), float tol,
        float *xmin);

float df;
float prob_thresh;
 
main(int argc, char *argv[])
{
    float t;
    float min;

/* ---------------------------------------------------------------------------- */
/*  CHECK TO SEE IF CORRECT NUMBER OF ARGUMENTS WERE PASSED.                    */

    if(argc<3) {
        printf("t: need probability threshold and degrees of freedom.\n");
        printf("   Will return TWO-sided t-test threshold.\n");
        exit(1);
    }

/* ---------------------------------------------------------------------------- */
/*  READ IN COMMAND LINE ARGUMENTS.                                             */

    prob_thresh = atof(argv[1]);
    df     = atof(argv[2]);
/*    printf("prob_thresh = %f\n",prob_thresh);
    printf("df          = %4.0f\n",df); */

    min=golden(0,3,50,(float (*) (float)) t2,0.0001,&t);

/*    printf("For p=%f, df=%d, threshold 2-tailed t=%f\n",prob_thresh,(int)df,t);
    printf("Deviation = %f\n",min); */
    printf("%8.5f\n",t);
}

/* ---------------------------------------------------------------------------- */
/* t distribution.                                                              */
float t2(float t_test)
{
/*    printf("prob_thresh = %f, current t_prob = %f, t = %f\n",prob_thresh,
        1.-betai(df/2.,0.5,(df/(df+(t_test*t_test)))),t_test);
    printf("%f\n",(float) fabs((double) (prob_thresh-betai(df/2.,0.5,(df/(df+(t_test*t_test))))))); */
    return (float) fabs(prob_thresh-betai(df/2.,0.5,(df/(df+(t_test*t_test)))));
}
