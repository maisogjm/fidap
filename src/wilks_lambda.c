#include <math.h>
#include <stdio.h>
#include <stdlib.h>

float lambda(float wilks_lambda);
float betai(float a, float b, float x);
float golden(float ax, float bx, float cx, float (*f)(float), float tol,
        float *xmin);
float df1,df2;                /* Degrees of freedom for the F-test.             */
float prob_thresh;            /* Probabilidy threshold.                         */
float wilks_lambda_power;     /* Wilks' Lambda raised to the (1/t)-th power.    */
float w,t;                    /* Intermediate variables for converting to F-test*/
float df_ratio;               /* df2/df1                                        */

main(int argc, char *argv[])
{
    /* 01.20.96 : JM Maisog
        Converts a Wilks' Lambda map to an F-test map, using an approximate
        transformation.  See Rencher AC, Methods of Multivariate Analysis,
        New York: John Wiley & Sons, 1995, eq. (6.13), p. 182.                  */

    int num_voxX4;            /* Number of bytes in image.                      */
    int num_vox;              /* Number of voxels in image.                     */
    int vox_index;            /* Index into voxels in image.                    */

    int read_status;          /* Check on read status of fread.                 */
    float wilks_lambda;       /* Wilks' Lambda.                                 */
    float p, nu_H, nu_E;      /* Parameters of Wilks' Lambda.                   */

    float min;


/* ---------------------------------------------------------------------------- */
/*  CHECK TO SEE IF CORRECT NUMBER OF ARGUMENTS WERE PASSED.                    */

    if (argc!=5) {
        printf("wilks_lambda: need significance threshold, p, nu-H, and nu-E.\n\n");
        printf("Will return Wilk's lambda for the given significance threshold,\n");
        printf("using an approximate Wilks' Lambda-to-F-test transform.\n");
        printf("This program may hang if nu-E is too big.  For p=1, consider using\n");
        printf("wilks_lambda_p1, which is an exact Wilks' Lambda-to-F-test transform.\n");

        exit(1);
    }

/* ---------------------------------------------------------------------------- */
/*  READ IN COMMAND LINE ARGUMENTS.                                             */

    prob_thresh=atof(argv[1]);
    p    = atof(argv[2]);
    nu_H = atof(argv[3]);
    nu_E = atof(argv[4]);
/*
    printf("p    = %f\n",p);
    printf("nu_H = %f\n",nu_H);
    printf("nu_E = %f\n",nu_E);
 */


/* ---------------------------------------------------------------------------- */
/* COMPUTE df1, df2, w, AND t.                                                  */

    t        = sqrt( ((p*p*nu_H*nu_H)-4.) / ((p*p)+(nu_H*nu_H)-5.) );
    w        = nu_E+nu_H-(0.5*(p+nu_H+1.));
    df1      = p*nu_H;
    df2      = (w*t)-(0.5*((p*nu_H)-2.));
    df_ratio = df2/df1;
/*
    printf("t    = %f\n",t);
    printf("w    = %f\n",w);
    printf("df1  = %f\n",df1);
    printf("df2  = %f\n",df2);
*/

/* ---------------------------------------------------------------------------- */
/* CONVERT WILKS' LAMBDA TO APPROXIMATE F-TEST.                                 */
    min=golden(0.,0.5,1.,(float (*) (float)) lambda,0.0001,&wilks_lambda);
    printf("%g\n",wilks_lambda);
}

/* ---------------------------------------------------------------------------- */
/* Approximate Wilks'Lambda distribution.                                       */
float lambda(float wilks_lambda)
{
    float F_test;

    wilks_lambda_power=pow(wilks_lambda,1./t);
    F_test=((1.-wilks_lambda_power)/wilks_lambda_power)*df_ratio;

/*
    printf("wilks_lambda=%g, wilks_lambda_power=%g, F=%g\n",wilks_lambda,wilks_lambda_power,F_test);
*/
    return (float) fabs(prob_thresh-betai(df2/2.,df1/2,(df2/(df2+(df1*F_test)))));

}

